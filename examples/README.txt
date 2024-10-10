
## Sentiment Analysis with SVM

This project implements a sentiment analysis model using a Support Vector Machine (SVM) in Python. The model is trained on the **20 Newsgroups Dataset** using `scikit-learn`, with text preprocessing handled by the `TfidfVectorizer`.

### Table of Contents

1. [Installation](#installation)
2. [Running the Sentiment Analysis](#running-the-sentiment-analysis)
3. [Project Structure](#project-structure)
4. [How the Code Works](#how-the-code-works)
5. [Requirements](#requirements)

---

### Installation

To get started, clone the repository and set up the required dependencies. Ensure you are working in a virtual environment:

```bash
git clone https://github.com/YOUR_USERNAME/scikit-learn.git
cd scikit-learn
python3 -m venv env
source env/bin/activate  # For Mac/Linux
env\Scriptsctivate      # For Windows
pip install -r requirements.txt  # If you have one, otherwise install manually
```

Install dependencies manually if no `requirements.txt` exists:
```bash
pip install numpy scipy scikit-learn matplotlib pandas
```

### Running the Sentiment Analysis

1. After the installation, run the sentiment analysis script located in the `examples` directory:
   ```bash
   python examples/sentiment_analysis.py
   ```

2. The code will:
   - Load the **20 Newsgroups** dataset for specific categories.
   - Preprocess the text data using `TfidfVectorizer`.
   - Train the SVM model on the training data.
   - Test the model and display classification metrics.

### Project Structure

```
scikit-learn/
│
├── examples/
│   └── sentiment_analysis.py   # The main script for SVM-based sentiment analysis
├── README.md                   # Project documentation (this file)
├── env/                        # Virtual environment (after setup)
└── ...
```

### How the Code Works

1. **Data Loading**: The dataset used in the project is the `20 Newsgroups` dataset provided by `scikit-learn`.
2. **Text Preprocessing**: The text data is transformed into TF-IDF features using `TfidfVectorizer`.
3. **Model Training**: A Support Vector Machine (`SVC`) is used to classify the text data.
4. **Evaluation**: The performance of the model is evaluated using a test set, and the results are displayed using `metrics.classification_report`.

### Requirements

- Python 3.x
- numpy
- scipy
- scikit-learn
- matplotlib
- pandas
