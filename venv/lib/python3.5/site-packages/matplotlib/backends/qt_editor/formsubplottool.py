from matplotlib.backends.qt_compat import QtWidgets


class UiSubplotTool(QtWidgets.QDialog):

    def __init__(self, *args, **kwargs):
        super(UiSubplotTool, self).__init__(*args, **kwargs)
        self.setObjectName("SubplotTool")
        self._widgets = {}

        layout = QtWidgets.QHBoxLayout()
        self.setLayout(layout)

        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left)
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right)

        box = QtWidgets.QGroupBox("Borders")
        left.addWidget(box)
        inner = QtWidgets.QFormLayout(box)
        for side in ["top", "bottom", "left", "right"]:
            self._widgets[side] = widget = QtWidgets.QDoubleSpinBox()
            widget.setMinimum(0)
            widget.setMaximum(1)
            widget.setDecimals(3)
            widget.setSingleStep(.005)
            widget.setKeyboardTracking(False)
            inner.addRow(side, widget)
        left.addStretch(1)

        box = QtWidgets.QGroupBox("Spacings")
        right.addWidget(box)
        inner = QtWidgets.QFormLayout(box)
        for side in ["hspace", "wspace"]:
            self._widgets[side] = widget = QtWidgets.QDoubleSpinBox()
            widget.setMinimum(0)
            widget.setMaximum(1)
            widget.setDecimals(3)
            widget.setSingleStep(.005)
            widget.setKeyboardTracking(False)
            inner.addRow(side, widget)
        right.addStretch(1)

        widget = QtWidgets.QPushButton("Export values")
        self._widgets["Export values"] = widget
        # Don't trigger on <enter>, which is used to input values.
        widget.setAutoDefault(False)
        left.addWidget(widget)

        for action in ["Tight layout", "Reset", "Close"]:
            self._widgets[action] = widget = QtWidgets.QPushButton(action)
            widget.setAutoDefault(False)
            right.addWidget(widget)

        self._widgets["Close"].setFocus()
