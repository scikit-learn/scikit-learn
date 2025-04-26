import tracemalloc
import os
import pytest
import numpy as np
import re
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LinearRegression
from sklearn.metrics import confusion_matrix
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics._plot.confusion_matrix_sankey import ConfusionMatrixSankeyDisplay

def test_initialization():
    # 测试初始化是否正确设置属性
    cm = np.array([[2, 1], [0, 3]])
    display_labels = ['A', 'B']
    save_path = "test.html"
    sankey = ConfusionMatrixSankeyDisplay(
        cm, display_labels=display_labels, save_path=save_path
    )
    assert np.array_equal(sankey.cm, cm)
    assert sankey.display_labels == display_labels
    assert sankey.save_path == save_path

def test_from_predictions():
    # 测试 from_predictions 生成的混淆矩阵是否正确
    y_true = [0, 1, 0, 1]
    y_pred = [0, 1, 1, 0]
    display = ConfusionMatrixSankeyDisplay.from_predictions(y_true, y_pred)
    expected_cm = confusion_matrix(y_true, y_pred)
    assert np.array_equal(display.cm, expected_cm)

def test_from_predictions_normalize():
    # 测试归一化后的混淆矩阵是否正确
    y_true = [0, 0, 1, 1]
    y_pred = [0, 0, 1, 1]
    display = ConfusionMatrixSankeyDisplay.from_predictions(
        y_true, y_pred, normalize='true'
    )
    raw_cm = confusion_matrix(y_true, y_pred, normalize='true')
    expected_cm = (raw_cm * 100).astype(int)
    assert np.array_equal(display.cm, expected_cm)

def test_from_estimator():
    # 测试 from_estimator 是否正确使用分类器的预测结果
    X, y = load_iris(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    clf = DummyClassifier(strategy='most_frequent')
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    display = ConfusionMatrixSankeyDisplay.from_estimator(clf, X_test, y_test)
    expected_cm = confusion_matrix(y_test, y_pred)
    assert np.array_equal(display.cm, expected_cm)

def test_from_estimator_non_classifier():
    # 测试传入非分类器时是否抛出错误
    reg = LinearRegression()
    X = [[0], [1], [2], [3]]
    y = [0, 1, 0, 1]
    reg.fit(X, y)
    with pytest.raises(ValueError):
        ConfusionMatrixSankeyDisplay.from_estimator(reg, X, y)

def test_plot_saves_file(tmpdir):
    # 测试 plot 方法是否生成文件
    save_path = tmpdir.join("sankey.html")
    cm = np.array([[2, 0], [0, 3]])
    display = ConfusionMatrixSankeyDisplay(cm, save_path=str(save_path))
    display.plot()
    assert save_path.exists()

def test_hex_to_rgb():
    # 测试颜色转换方法是否正确
    sankey = ConfusionMatrixSankeyDisplay([[1]], display_labels=['A'])
    assert sankey._hex_to_rgb("#ff0000") == "255,0,0"
    assert sankey._hex_to_rgb("#00ff00") == "0,255,0"
    assert sankey._hex_to_rgb("#0000ff") == "0,0,255"

def test_display_labels_default():
    # 测试未提供 display_labels 时是否使用默认值
    cm = np.array([[1, 0], [0, 1]])
    sankey = ConfusionMatrixSankeyDisplay(cm)
    assert np.array_equal(sankey.display_labels, [0, 1])

def _extract_title(text):
    match = re.search(r'"title":\s*{\s*"text":\s*"([^"]+)"', text)
    return match.group(1) if match else ""

def _parse_html_content(html_path):
    """Parse HTML content for validation"""
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    return {
        "has_sankey": "Sankey" in content,
        "title": _extract_title(content),
        "node_count": content.count('"node":{'),
        "link_count": content.count('"link":{'),
        "cm_values": _extract_cm_values(content),
        "hover_template": _extract_hover_template(content),
        "link_colors": _extract_link_colors(content)
    }


def _extract_between(text, start, end):
    parts = text.split(start)
    return parts[1].split(end)[0] if len(parts) > 1 else ""


def _extract_cm_values(html):
    link_value_match = re.search(
        r'"link":\s*\{[^}]+"value":\s*\[([\d\s,]+)\]',
        html,
        re.DOTALL
    )
    if link_value_match:
        values_str = link_value_match.group(1)
        values = [
            int(float(v.strip()))
            for v in re.split(r'\s*,\s*', values_str.replace('\n', ''))
            if v.strip()
        ]
        return values
    return []


def _extract_hover_template(html):
    parts = html.split('"hovertemplate":')
    return parts[1].split('",')[0] if len(parts) > 1 else ""


def _extract_link_colors(html):
    return re.findall(r'"color": "(rgba?\([^)]+\))",', html)


# Test Cases

def test_html_content(tmpdir):
    save_path = tmpdir.join("sankey.html")
    cm = np.array([[2, 1], [3, 4]])
    display = ConfusionMatrixSankeyDisplay(
        cm,
        display_labels=["Cat", "Dog"],
        save_path=str(save_path)
    ).plot()

    content = _parse_html_content(save_path)

    assert content["has_sankey"]
    assert "Confusion Matrix Sankey Diagram" in content["title"]
    assert content["node_count"] == 1
    assert content["link_count"] == 1


def test_normalized_html_content(tmpdir):
    save_path = tmpdir.join("sankey_norm.html")
    cm = np.array([[50, 50], [30, 70]])
    display = ConfusionMatrixSankeyDisplay(
        cm,
        display_labels=["A", "B"],
        save_path=str(save_path)
    ).plot(values_format=".1f%")

    content = _parse_html_content(save_path)
    assert "Percentage: %{value:.1f}%" in content["hover_template"]


@pytest.mark.benchmark
def test_large_scale_performance(tmpdir, benchmark):
    """Benchmark testing for large-scale confusion matrix"""
    # 生成1000x1000的随机混淆矩阵（稀疏）
    np.random.seed(42)
    n_classes = 500
    cm = np.random.randint(0, 10, size=(n_classes, n_classes))

    # 创建显示对象
    display = ConfusionMatrixSankeyDisplay(confusion_matrix=cm,
        display_labels=np.arange(n_classes),
        save_path=str(tmpdir.join("large_sankey.html")))

    # 定义性能测试函数

    def _run():
        # 内存使用监控
        tracemalloc.start()
        display.plot(width=1600, height=1200)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return peak / 1024 ** 2  # 转换为MB

    # 执行基准测试
    peak_mem = benchmark.pedantic(_run, rounds=3, warmup_rounds=1)

    # 断言内存峰值不超过500MB（根据实际硬件调整阈值）
    assert peak_mem < 500, f"内存使用过高：{peak_mem:.2f}MB"

    # 验证生成的文件
    assert os.path.exists(display.save_path)
    file_size = os.path.getsize(display.save_path) / 1024 ** 2
    assert file_size < 50, f"生成文件过大：{file_size:.2f}MB"

