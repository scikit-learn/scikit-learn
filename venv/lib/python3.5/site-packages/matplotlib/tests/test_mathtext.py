from __future__ import absolute_import, division, print_function

import six

import io
import re

import numpy as np
import pytest

import matplotlib
from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt
from matplotlib import mathtext


math_tests = [
    r'$a+b+\dot s+\dot{s}+\ldots$',
    r'$x \doteq y$',
    r'\$100.00 $\alpha \_$',
    r'$\frac{\$100.00}{y}$',
    r'$x   y$',
    r'$x+y\ x=y\ x<y\ x:y\ x,y\ x@y$',
    r'$100\%y\ x*y\ x/y x\$y$',
    r'$x\leftarrow y\ x\forall y\ x-y$',
    r'$x \sf x \bf x {\cal X} \rm x$',
    r'$x\ x\,x\;x\quad x\qquad x\!x\hspace{ 0.5 }y$',
    r'$\{ \rm braces \}$',
    r'$\left[\left\lfloor\frac{5}{\frac{\left(3\right)}{4}} y\right)\right]$',
    r'$\left(x\right)$',
    r'$\sin(x)$',
    r'$x_2$',
    r'$x^2$',
    r'$x^2_y$',
    r'$x_y^2$',
    r'$\prod_{i=\alpha_{i+1}}^\infty$',
    r'$x = \frac{x+\frac{5}{2}}{\frac{y+3}{8}}$',
    r'$dz/dt = \gamma x^2 + {\rm sin}(2\pi y+\phi)$',
    r'Foo: $\alpha_{i+1}^j = {\rm sin}(2\pi f_j t_i) e^{-5 t_i/\tau}$',
    r'$\mathcal{R}\prod_{i=\alpha_{i+1}}^\infty a_i \sin(2 \pi f x_i)$',
    r'Variable $i$ is good',
    r'$\Delta_i^j$',
    r'$\Delta^j_{i+1}$',
    r'$\ddot{o}\acute{e}\grave{e}\hat{O}\breve{\imath}\tilde{n}\vec{q}$',
    r"$\arccos((x^i))$",
    r"$\gamma = \frac{x=\frac{6}{8}}{y} \delta$",
    r'$\limsup_{x\to\infty}$',
    r'$\oint^\infty_0$',
    r"$f'\quad f'''(x)\quad ''/\mathrm{yr}$",
    r'$\frac{x_2888}{y}$',
    r"$\sqrt[3]{\frac{X_2}{Y}}=5$",
    r"$\sqrt[5]{\prod^\frac{x}{2\pi^2}_\infty}$",
    r"$\sqrt[3]{x}=5$",
    r'$\frac{X}{\frac{X}{Y}}$',
    r"$W^{3\beta}_{\delta_1 \rho_1 \sigma_2} = U^{3\beta}_{\delta_1 \rho_1} + \frac{1}{8 \pi 2} \int^{\alpha_2}_{\alpha_2} d \alpha^\prime_2 \left[\frac{ U^{2\beta}_{\delta_1 \rho_1} - \alpha^\prime_2U^{1\beta}_{\rho_1 \sigma_2} }{U^{0\beta}_{\rho_1 \sigma_2}}\right]$",
    r'$\mathcal{H} = \int d \tau \left(\epsilon E^2 + \mu H^2\right)$',
    r'$\widehat{abc}\widetilde{def}$',
    '$\\Gamma \\Delta \\Theta \\Lambda \\Xi \\Pi \\Sigma \\Upsilon \\Phi \\Psi \\Omega$',
    '$\\alpha \\beta \\gamma \\delta \\epsilon \\zeta \\eta \\theta \\iota \\lambda \\mu \\nu \\xi \\pi \\kappa \\rho \\sigma \\tau \\upsilon \\phi \\chi \\psi$',

    # The examples prefixed by 'mmltt' are from the MathML torture test here:
        # http://www.mozilla.org/projects/mathml/demo/texvsmml.xhtml
    r'${x}^{2}{y}^{2}$',
    r'${}_{2}F_{3}$',
    r'$\frac{x+{y}^{2}}{k+1}$',
    r'$x+{y}^{\frac{2}{k+1}}$',
    r'$\frac{a}{b/2}$',
    r'${a}_{0}+\frac{1}{{a}_{1}+\frac{1}{{a}_{2}+\frac{1}{{a}_{3}+\frac{1}{{a}_{4}}}}}$',
    r'${a}_{0}+\frac{1}{{a}_{1}+\frac{1}{{a}_{2}+\frac{1}{{a}_{3}+\frac{1}{{a}_{4}}}}}$',
    r'$\binom{n}{k/2}$',
    r'$\binom{p}{2}{x}^{2}{y}^{p-2}-\frac{1}{1-x}\frac{1}{1-{x}^{2}}$',
    r'${x}^{2y}$',
    r'$\sum _{i=1}^{p}\sum _{j=1}^{q}\sum _{k=1}^{r}{a}_{ij}{b}_{jk}{c}_{ki}$',
    r'$\sqrt{1+\sqrt{1+\sqrt{1+\sqrt{1+\sqrt{1+\sqrt{1+\sqrt{1+x}}}}}}}$',
    r'$\left(\frac{{\partial }^{2}}{\partial {x}^{2}}+\frac{{\partial }^{2}}{\partial {y}^{2}}\right){|\varphi \left(x+iy\right)|}^{2}=0$',
    r'${2}^{{2}^{{2}^{x}}}$',
    r'${\int }_{1}^{x}\frac{\mathrm{dt}}{t}$',
    r'$\int {\int }_{D}\mathrm{dx} \mathrm{dy}$',
    # mathtex doesn't support array
    # 'mmltt18'    : r'$f\left(x\right)=\left\{\begin{array}{cc}\hfill 1/3\hfill & \text{if_}0\le x\le 1;\hfill \\ \hfill 2/3\hfill & \hfill \text{if_}3\le x\le 4;\hfill \\ \hfill 0\hfill & \text{elsewhere.}\hfill \end{array}$',
    # mathtex doesn't support stackrel
    # 'mmltt19'    : ur'$\stackrel{\stackrel{k\text{times}}{\ufe37}}{x+...+x}$',
    r'${y}_{{x}^{2}}$',
    # mathtex doesn't support the "\text" command
    # 'mmltt21'    : r'$\sum _{p\text{\prime}}f\left(p\right)={\int }_{t>1}f\left(t\right) d\pi \left(t\right)$',
    # mathtex doesn't support array
    # 'mmltt23'    : r'$\left(\begin{array}{cc}\hfill \left(\begin{array}{cc}\hfill a\hfill & \hfill b\hfill \\ \hfill c\hfill & \hfill d\hfill \end{array}\right)\hfill & \hfill \left(\begin{array}{cc}\hfill e\hfill & \hfill f\hfill \\ \hfill g\hfill & \hfill h\hfill \end{array}\right)\hfill \\ \hfill 0\hfill & \hfill \left(\begin{array}{cc}\hfill i\hfill & \hfill j\hfill \\ \hfill k\hfill & \hfill l\hfill \end{array}\right)\hfill \end{array}\right)$',
    # mathtex doesn't support array
    # 'mmltt24'   : u'$det|\\begin{array}{ccccc}\\hfill {c}_{0}\\hfill & \\hfill {c}_{1}\\hfill & \\hfill {c}_{2}\\hfill & \\hfill \\dots \\hfill & \\hfill {c}_{n}\\hfill \\\\ \\hfill {c}_{1}\\hfill & \\hfill {c}_{2}\\hfill & \\hfill {c}_{3}\\hfill & \\hfill \\dots \\hfill & \\hfill {c}_{n+1}\\hfill \\\\ \\hfill {c}_{2}\\hfill & \\hfill {c}_{3}\\hfill & \\hfill {c}_{4}\\hfill & \\hfill \\dots \\hfill & \\hfill {c}_{n+2}\\hfill \\\\ \\hfill \\u22ee\\hfill & \\hfill \\u22ee\\hfill & \\hfill \\u22ee\\hfill & \\hfill \\hfill & \\hfill \\u22ee\\hfill \\\\ \\hfill {c}_{n}\\hfill & \\hfill {c}_{n+1}\\hfill & \\hfill {c}_{n+2}\\hfill & \\hfill \\dots \\hfill & \\hfill {c}_{2n}\\hfill \\end{array}|>0$',
    r'${y}_{{x}_{2}}$',
    r'${x}_{92}^{31415}+\pi $',
    r'${x}_{{y}_{b}^{a}}^{{z}_{c}^{d}}$',
    r'${y}_{3}^{\prime \prime \prime }$',
    r"$\left( \xi \left( 1 - \xi \right) \right)$", # Bug 2969451
    r"$\left(2 \, a=b\right)$", # Sage bug #8125
    r"$? ! &$", # github issue #466
    r'$\operatorname{cos} x$', # github issue #553
    r'$\sum _{\genfrac{}{}{0}{}{0\leq i\leq m}{0<j<n}}P\left(i,j\right)$',
    r"$\left\Vert a \right\Vert \left\vert b \right\vert \left| a \right| \left\| b\right\| \Vert a \Vert \vert b \vert$",
    r'$\mathring{A}  \stackrel{\circ}{A}  \AA$',
    r'$M \, M \thinspace M \/ M \> M \: M \; M \ M \enspace M \quad M \qquad M \! M$',
    r'$\Cup$ $\Cap$ $\leftharpoonup$ $\barwedge$ $\rightharpoonup$',
    r'$\dotplus$ $\doteq$ $\doteqdot$ $\ddots$',
    r'$xyz^kx_kx^py^{p-2} d_i^jb_jc_kd x^j_i E^0 E^0_u$', # github issue #4873
    r'${xyz}^k{x}_{k}{x}^{p}{y}^{p-2} {d}_{i}^{j}{b}_{j}{c}_{k}{d} {x}^{j}_{i}{E}^{0}{E}^0_u$',
    r'${\int}_x^x x\oint_x^x x\int_{X}^{X}x\int_x x \int^x x \int_{x} x\int^{x}{\int}_{x} x{\int}^{x}_{x}x$',
    r'testing$^{123}$',
    ' '.join('$\\' + p + '$' for p in sorted(mathtext.Parser._snowflake)),
    r'$6-2$; $-2$; $ -2$; ${-2}$; ${  -2}$; $20^{+3}_{-2}$',
    r'$\overline{\omega}^x \frac{1}{2}_0^x$', # github issue #5444
    r'$,$ $.$ $1{,}234{, }567{ , }890$ and $1,234,567,890$', # github issue 5799
    r'$\left(X\right)_{a}^{b}$', # github issue 7615
    r'$\dfrac{\$100.00}{y}$', # github issue #1888
]

digits = "0123456789"
uppercase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
lowercase = "abcdefghijklmnopqrstuvwxyz"
uppergreek = ("\\Gamma \\Delta \\Theta \\Lambda \\Xi \\Pi \\Sigma \\Upsilon \\Phi \\Psi "
              "\\Omega")
lowergreek = ("\\alpha \\beta \\gamma \\delta \\epsilon \\zeta \\eta \\theta \\iota "
              "\\lambda \\mu \\nu \\xi \\pi \\kappa \\rho \\sigma \\tau \\upsilon "
              "\\phi \\chi \\psi")
all = [digits, uppercase, lowercase, uppergreek, lowergreek]

font_test_specs = [
    ([], all),
    (['mathrm'], all),
    (['mathbf'], all),
    (['mathit'], all),
    (['mathtt'], [digits, uppercase, lowercase]),
    (['mathcircled'], [digits, uppercase, lowercase]),
    (['mathrm', 'mathcircled'], [digits, uppercase, lowercase]),
    (['mathbf', 'mathcircled'], [digits, uppercase, lowercase]),
    (['mathbb'], [digits, uppercase, lowercase,
                  r'\Gamma \Pi \Sigma \gamma \pi']),
    (['mathrm', 'mathbb'], [digits, uppercase, lowercase,
                            r'\Gamma \Pi \Sigma \gamma \pi']),
    (['mathbf', 'mathbb'], [digits, uppercase, lowercase,
                            r'\Gamma \Pi \Sigma \gamma \pi']),
    (['mathcal'], [uppercase]),
    (['mathfrak'], [uppercase, lowercase]),
    (['mathbf', 'mathfrak'], [uppercase, lowercase]),
    (['mathscr'], [uppercase, lowercase]),
    (['mathsf'], [digits, uppercase, lowercase]),
    (['mathrm', 'mathsf'], [digits, uppercase, lowercase]),
    (['mathbf', 'mathsf'], [digits, uppercase, lowercase])
    ]

font_tests = []
for fonts, chars in font_test_specs:
    wrapper = [' '.join(fonts), ' $']
    for font in fonts:
        wrapper.append(r'\%s{' % font)
    wrapper.append('%s')
    for font in fonts:
        wrapper.append('}')
    wrapper.append('$')
    wrapper = ''.join(wrapper)

    for set in chars:
        font_tests.append(wrapper % set)


@pytest.fixture
def baseline_images(request, fontset, index):
    return ['%s_%s_%02d' % (request.param, fontset, index)]


@pytest.mark.parametrize('index, test', enumerate(math_tests),
                         ids=[str(index) for index in range(len(math_tests))])
@pytest.mark.parametrize('fontset',
                         ['cm', 'stix', 'stixsans', 'dejavusans',
                          'dejavuserif'])
@pytest.mark.parametrize('baseline_images', ['mathtext'], indirect=True)
@image_comparison(baseline_images=None)
def test_mathtext_rendering(baseline_images, fontset, index, test):
    matplotlib.rcParams['mathtext.fontset'] = fontset
    fig = plt.figure(figsize=(5.25, 0.75))
    fig.text(0.5, 0.5, test,
             horizontalalignment='center', verticalalignment='center')


@pytest.mark.parametrize('index, test', enumerate(font_tests),
                         ids=[str(index) for index in range(len(font_tests))])
@pytest.mark.parametrize('fontset',
                         ['cm', 'stix', 'stixsans', 'dejavusans',
                          'dejavuserif'])
@pytest.mark.parametrize('baseline_images', ['mathfont'], indirect=True)
@image_comparison(baseline_images=None, extensions=['png'])
def test_mathfont_rendering(baseline_images, fontset, index, test):
    matplotlib.rcParams['mathtext.fontset'] = fontset
    fig = plt.figure(figsize=(5.25, 0.75))
    fig.text(0.5, 0.5, test,
             horizontalalignment='center', verticalalignment='center')


def test_fontinfo():
    import matplotlib.font_manager as font_manager
    import matplotlib.ft2font as ft2font
    fontpath = font_manager.findfont("DejaVu Sans")
    font = ft2font.FT2Font(fontpath)
    table = font.get_sfnt_table("head")
    assert table['version'] == (1, 0)


@pytest.mark.parametrize(
    'math, msg',
    [
        (r'$\hspace{}$', r'Expected \hspace{n}'),
        (r'$\hspace{foo}$', r'Expected \hspace{n}'),
        (r'$\frac$', r'Expected \frac{num}{den}'),
        (r'$\frac{}{}$', r'Expected \frac{num}{den}'),
        (r'$\stackrel$', r'Expected \stackrel{num}{den}'),
        (r'$\stackrel{}{}$', r'Expected \stackrel{num}{den}'),
        (r'$\binom$', r'Expected \binom{num}{den}'),
        (r'$\binom{}{}$', r'Expected \binom{num}{den}'),
        (r'$\genfrac$',
         r'Expected \genfrac{ldelim}{rdelim}{rulesize}{style}{num}{den}'),
        (r'$\genfrac{}{}{}{}{}{}$',
         r'Expected \genfrac{ldelim}{rdelim}{rulesize}{style}{num}{den}'),
        (r'$\sqrt$', r'Expected \sqrt{value}'),
        (r'$\sqrt f$', r'Expected \sqrt{value}'),
        (r'$\overline$', r'Expected \overline{value}'),
        (r'$\overline{}$', r'Expected \overline{value}'),
        (r'$\leftF$', r'Expected a delimiter'),
        (r'$\rightF$', r'Unknown symbol: \rightF'),
        (r'$\left(\right$', r'Expected a delimiter'),
        (r'$\left($', r'Expected "\right"'),
        (r'$\dfrac$', r'Expected \dfrac{num}{den}'),
        (r'$\dfrac{}{}$', r'Expected \dfrac{num}{den}'),
    ],
    ids=[
        'hspace without value',
        'hspace with invalid value',
        'frac without parameters',
        'frac with empty parameters',
        'stackrel without parameters',
        'stackrel with empty parameters',
        'binom without parameters',
        'binom with empty parameters',
        'genfrac without parameters',
        'genfrac with empty parameters',
        'sqrt without parameters',
        'sqrt with invalid value',
        'overline without parameters',
        'overline with empty parameter',
        'left with invalid delimiter',
        'right with invalid delimiter',
        'unclosed parentheses with sizing',
        'unclosed parentheses without sizing',
        'dfrac without parameters',
        'dfrac with empty parameters',
    ]
)
def test_mathtext_exceptions(math, msg):
    parser = mathtext.MathTextParser('agg')

    with pytest.raises(ValueError) as excinfo:
        parser.parse(math)
    excinfo.match(re.escape(msg))


def test_single_minus_sign():
    plt.figure(figsize=(0.3, 0.3))
    plt.text(0.5, 0.5, '$-$')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])

    buff = io.BytesIO()
    plt.savefig(buff, format="rgba", dpi=1000)
    array = np.fromstring(buff.getvalue(), dtype=np.uint8)

    # If this fails, it would be all white
    assert not np.all(array == 0xff)
