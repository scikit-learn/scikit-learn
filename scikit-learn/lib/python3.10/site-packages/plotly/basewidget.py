from copy import deepcopy
import pathlib
from traitlets import List, Dict, observe, Integer
from plotly.io._renderers import display_jupyter_version_warnings

from .basedatatypes import BaseFigure, BasePlotlyType
from .callbacks import BoxSelector, LassoSelector, InputDeviceState, Points
from .serializers import custom_serializers
import anywidget


class BaseFigureWidget(BaseFigure, anywidget.AnyWidget):
    """
    Base class for FigureWidget. The FigureWidget class is code-generated as a
    subclass
    """

    _esm = pathlib.Path(__file__).parent / "package_data" / "widgetbundle.js"

    # ### _data and _layout ###
    # These properties store the current state of the traces and
    # layout as JSON-style dicts. These dicts do not store any subclasses of
    # `BasePlotlyType`
    #
    # Note: These are only automatically synced with the frontend on full
    # assignment, not on mutation. We use this fact to only directly sync
    # them to the front-end on FigureWidget construction. All other updates
    # are made using mutation, and they are manually synced to the frontend
    # using the relayout/restyle/update/etc. messages.
    _widget_layout = Dict().tag(sync=True, **custom_serializers)
    _widget_data = List().tag(sync=True, **custom_serializers)
    _config = Dict().tag(sync=True, **custom_serializers)

    # ### Python -> JS message properties ###
    # These properties are used to send messages from Python to the
    # frontend. Messages are sent by assigning the message contents to the
    # appropriate _py2js_* property and then immediatly assigning None to the
    # property.
    #
    # See JSDoc comments in the FigureModel class in js/src/Figure.js for
    # detailed descriptions of the messages.
    _py2js_addTraces = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _py2js_restyle = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _py2js_relayout = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _py2js_update = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _py2js_animate = Dict(allow_none=True).tag(sync=True, **custom_serializers)

    _py2js_deleteTraces = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _py2js_moveTraces = Dict(allow_none=True).tag(sync=True, **custom_serializers)

    _py2js_removeLayoutProps = Dict(allow_none=True).tag(
        sync=True, **custom_serializers
    )
    _py2js_removeTraceProps = Dict(allow_none=True).tag(sync=True, **custom_serializers)

    # ### JS -> Python message properties ###
    # These properties are used to receive messages from the frontend.
    # Messages are received by defining methods that observe changes to these
    # properties. Receive methods are named `_handler_js2py_*` where '*' is
    # the name of the corresponding message property.  Receive methods are
    # responsible for setting the message property to None after retreiving
    # the message data.
    #
    # See JSDoc comments in the FigureModel class in js/src/Figure.js for
    # detailed descriptions of the messages.
    _js2py_traceDeltas = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _js2py_layoutDelta = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _js2py_restyle = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _js2py_relayout = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _js2py_update = Dict(allow_none=True).tag(sync=True, **custom_serializers)
    _js2py_pointsCallback = Dict(allow_none=True).tag(sync=True, **custom_serializers)

    # ### Message tracking properties ###
    # The _last_layout_edit_id and _last_trace_edit_id properties are used
    # to keep track of the edit id of the message that most recently
    # requested an update to the Figures layout or traces respectively.
    #
    # We track this information because we don't want to update the Figure's
    # default layout/trace properties (_layout_defaults, _data_defaults)
    # while edits are in process. This can lead to inconsistent property
    # states.
    _last_layout_edit_id = Integer(0).tag(sync=True)
    _last_trace_edit_id = Integer(0).tag(sync=True)

    _set_trace_uid = True
    _allow_disable_validation = False

    # Constructor
    # -----------
    def __init__(
        self, data=None, layout=None, frames=None, skip_invalid=False, **kwargs
    ):
        # Call superclass constructors
        # ----------------------------
        # Note: We rename layout to layout_plotly because to deconflict it
        # with the `layout` constructor parameter of the `widgets.DOMWidget`
        # ipywidgets class
        super(BaseFigureWidget, self).__init__(
            data=data,
            layout_plotly=layout,
            frames=frames,
            skip_invalid=skip_invalid,
            **kwargs,
        )

        # Validate Frames
        # ---------------
        # Frames are not supported by figure widget
        if self._frame_objs:
            BaseFigureWidget._display_frames_error()

        # Message States
        # --------------
        # ### Layout ###

        # _last_layout_edit_id is described above
        self._last_layout_edit_id = 0

        # _layout_edit_in_process is set to True if there are layout edit
        # operations that have been sent to the frontend that haven't
        # completed yet.
        self._layout_edit_in_process = False

        # _waiting_edit_callbacks is a list of callback functions that
        # should be executed as soon as all pending edit operations are
        # completed
        self._waiting_edit_callbacks = []

        # ### Trace ###
        # _last_trace_edit_id: described above
        self._last_trace_edit_id = 0

        # _trace_edit_in_process is set to True if there are trace edit
        # operations that have been sent to the frontend that haven't
        # completed yet.
        self._trace_edit_in_process = False

        # View count
        # ----------
        # ipywidget property that stores the number of active frontend
        # views of this widget
        self._view_count = 0

        # Initialize widget layout and data for third-party widget integration
        # --------------------------------------------------------------------
        self._widget_layout = deepcopy(self._layout_obj._props)
        self._widget_data = deepcopy(self._data)

    def show(self, *args, **kwargs):
        return self

    # Python -> JavaScript Messages
    # -----------------------------
    def _send_relayout_msg(self, layout_data, source_view_id=None):
        """
        Send Plotly.relayout message to the frontend

        Parameters
        ----------
        layout_data : dict
            Plotly.relayout layout data
        source_view_id : str
            UID of view that triggered this relayout operation
            (e.g. By the user clicking 'zoom' in the toolbar). None if the
            operation was not triggered by a frontend view
        """
        # Increment layout edit messages IDs
        # ----------------------------------
        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        # Build message
        # -------------
        msg_data = {
            "relayout_data": layout_data,
            "layout_edit_id": layout_edit_id,
            "source_view_id": source_view_id,
        }

        # Send message
        # ------------
        self._py2js_relayout = msg_data
        self._py2js_relayout = None

    def _send_restyle_msg(self, restyle_data, trace_indexes=None, source_view_id=None):
        """
        Send Plotly.restyle message to the frontend

        Parameters
        ----------
        restyle_data : dict
            Plotly.restyle restyle data
        trace_indexes : list[int]
            List of trace indexes that the restyle operation
            applies to
        source_view_id : str
            UID of view that triggered this restyle operation
            (e.g. By the user clicking the legend to hide a trace).
            None if the operation was not triggered by a frontend view
        """

        # Validate / normalize inputs
        # ---------------------------
        trace_indexes = self._normalize_trace_indexes(trace_indexes)

        # Increment layout/trace edit message IDs
        # ---------------------------------------
        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        trace_edit_id = self._last_trace_edit_id + 1
        self._last_trace_edit_id = trace_edit_id
        self._trace_edit_in_process = True

        # Build message
        # -------------
        restyle_msg = {
            "restyle_data": restyle_data,
            "restyle_traces": trace_indexes,
            "trace_edit_id": trace_edit_id,
            "layout_edit_id": layout_edit_id,
            "source_view_id": source_view_id,
        }

        # Send message
        # ------------
        self._py2js_restyle = restyle_msg
        self._py2js_restyle = None

    def _send_addTraces_msg(self, new_traces_data):
        """
        Send Plotly.addTraces message to the frontend

        Parameters
        ----------
        new_traces_data : list[dict]
            List of trace data for new traces as accepted by Plotly.addTraces
        """

        # Increment layout/trace edit message IDs
        # ---------------------------------------
        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        trace_edit_id = self._last_trace_edit_id + 1
        self._last_trace_edit_id = trace_edit_id
        self._trace_edit_in_process = True

        # Build message
        # -------------
        add_traces_msg = {
            "trace_data": new_traces_data,
            "trace_edit_id": trace_edit_id,
            "layout_edit_id": layout_edit_id,
        }

        # Send message
        # ------------
        self._py2js_addTraces = add_traces_msg
        self._py2js_addTraces = None

    def _send_moveTraces_msg(self, current_inds, new_inds):
        """
        Send Plotly.moveTraces message to the frontend

        Parameters
        ----------
        current_inds : list[int]
            List of current trace indexes
        new_inds : list[int]
            List of new trace indexes
        """

        # Build message
        # -------------
        move_msg = {"current_trace_inds": current_inds, "new_trace_inds": new_inds}

        # Send message
        # ------------
        self._py2js_moveTraces = move_msg
        self._py2js_moveTraces = None

    def _send_update_msg(
        self, restyle_data, relayout_data, trace_indexes=None, source_view_id=None
    ):
        """
        Send Plotly.update message to the frontend

        Parameters
        ----------
        restyle_data : dict
            Plotly.update restyle data
        relayout_data : dict
            Plotly.update relayout data
        trace_indexes : list[int]
            List of trace indexes that the update operation applies to
        source_view_id : str
            UID of view that triggered this update operation
            (e.g. By the user clicking a button).
            None if the operation was not triggered by a frontend view
        """

        # Validate / normalize inputs
        # ---------------------------
        trace_indexes = self._normalize_trace_indexes(trace_indexes)

        # Increment layout/trace edit message IDs
        # ---------------------------------------
        trace_edit_id = self._last_trace_edit_id + 1
        self._last_trace_edit_id = trace_edit_id
        self._trace_edit_in_process = True

        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        # Build message
        # -------------
        update_msg = {
            "style_data": restyle_data,
            "layout_data": relayout_data,
            "style_traces": trace_indexes,
            "trace_edit_id": trace_edit_id,
            "layout_edit_id": layout_edit_id,
            "source_view_id": source_view_id,
        }

        # Send message
        # ------------
        self._py2js_update = update_msg
        self._py2js_update = None

    def _send_animate_msg(
        self, styles_data, relayout_data, trace_indexes, animation_opts
    ):
        """
        Send Plotly.update message to the frontend

        Note: there is no source_view_id parameter because animations
        triggered by the fontend are not currently supported

        Parameters
        ----------
        styles_data : list[dict]
            Plotly.animate styles data
        relayout_data : dict
            Plotly.animate relayout data
        trace_indexes : list[int]
            List of trace indexes that the animate operation applies to
        """

        # Validate / normalize inputs
        # ---------------------------
        trace_indexes = self._normalize_trace_indexes(trace_indexes)

        # Increment layout/trace edit message IDs
        # ---------------------------------------
        trace_edit_id = self._last_trace_edit_id + 1
        self._last_trace_edit_id = trace_edit_id
        self._trace_edit_in_process = True

        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        # Build message
        # -------------
        animate_msg = {
            "style_data": styles_data,
            "layout_data": relayout_data,
            "style_traces": trace_indexes,
            "animation_opts": animation_opts,
            "trace_edit_id": trace_edit_id,
            "layout_edit_id": layout_edit_id,
            "source_view_id": None,
        }

        # Send message
        # ------------
        self._py2js_animate = animate_msg
        self._py2js_animate = None

    def _send_deleteTraces_msg(self, delete_inds):
        """
        Send Plotly.deleteTraces message to the frontend

        Parameters
        ----------
        delete_inds : list[int]
            List of trace indexes of traces to delete
        """

        # Increment layout/trace edit message IDs
        # ---------------------------------------
        trace_edit_id = self._last_trace_edit_id + 1
        self._last_trace_edit_id = trace_edit_id
        self._trace_edit_in_process = True

        layout_edit_id = self._last_layout_edit_id + 1
        self._last_layout_edit_id = layout_edit_id
        self._layout_edit_in_process = True

        # Build message
        # -------------
        delete_msg = {
            "delete_inds": delete_inds,
            "layout_edit_id": layout_edit_id,
            "trace_edit_id": trace_edit_id,
        }

        # Send message
        # ------------
        self._py2js_deleteTraces = delete_msg
        self._py2js_deleteTraces = None

    # JavaScript -> Python Messages
    # -----------------------------
    @observe("_js2py_traceDeltas")
    def _handler_js2py_traceDeltas(self, change):
        """
        Process trace deltas message from the frontend
        """

        # Receive message
        # ---------------
        msg_data = change["new"]
        if not msg_data:
            self._js2py_traceDeltas = None
            return

        trace_deltas = msg_data["trace_deltas"]
        trace_edit_id = msg_data["trace_edit_id"]

        # Apply deltas
        # ------------
        # We only apply the deltas if this message corresponds to the most
        # recent trace edit operation
        if trace_edit_id == self._last_trace_edit_id:
            # ### Loop over deltas ###
            for delta in trace_deltas:
                # #### Find existing trace for uid ###
                trace_uid = delta["uid"]
                trace_uids = [trace.uid for trace in self.data]
                trace_index = trace_uids.index(trace_uid)
                uid_trace = self.data[trace_index]

                # #### Transform defaults to delta ####
                delta_transform = BaseFigureWidget._transform_data(
                    uid_trace._prop_defaults, delta
                )

                # #### Remove overlapping properties ####
                # If a property is present in both _props and _prop_defaults
                # then we remove the copy from _props
                remove_props = self._remove_overlapping_props(
                    uid_trace._props, uid_trace._prop_defaults
                )

                # #### Notify frontend model of property removal ####
                if remove_props:
                    remove_trace_props_msg = {
                        "remove_trace": trace_index,
                        "remove_props": remove_props,
                    }
                    self._py2js_removeTraceProps = remove_trace_props_msg
                    self._py2js_removeTraceProps = None

                # #### Dispatch change callbacks ####
                self._dispatch_trace_change_callbacks(delta_transform, [trace_index])

            # ### Trace edits no longer in process ###
            self._trace_edit_in_process = False

            # ### Call any waiting trace edit callbacks ###
            if not self._layout_edit_in_process:
                while self._waiting_edit_callbacks:
                    self._waiting_edit_callbacks.pop()()

        self._js2py_traceDeltas = None

    @observe("_js2py_layoutDelta")
    def _handler_js2py_layoutDelta(self, change):
        """
        Process layout delta message from the frontend
        """

        # Receive message
        # ---------------
        msg_data = change["new"]
        if not msg_data:
            self._js2py_layoutDelta = None
            return

        layout_delta = msg_data["layout_delta"]
        layout_edit_id = msg_data["layout_edit_id"]

        # Apply delta
        # -----------
        # We only apply the delta if this message corresponds to the most
        # recent layout edit operation
        if layout_edit_id == self._last_layout_edit_id:
            # ### Transform defaults to delta ###
            delta_transform = BaseFigureWidget._transform_data(
                self._layout_defaults, layout_delta
            )

            # ### Remove overlapping properties ###
            # If a property is present in both _layout and _layout_defaults
            # then we remove the copy from _layout
            removed_props = self._remove_overlapping_props(
                self._widget_layout, self._layout_defaults
            )

            # ### Notify frontend model of property removal ###
            if removed_props:
                remove_props_msg = {"remove_props": removed_props}

                self._py2js_removeLayoutProps = remove_props_msg
                self._py2js_removeLayoutProps = None

            # ### Create axis objects ###
            # For example, when a SPLOM trace is created the layout defaults
            # may include axes that weren't explicitly defined by the user.
            for proppath in delta_transform:
                prop = proppath[0]
                match = self.layout._subplot_re_match(prop)
                if match and prop not in self.layout:
                    # We need to create a subplotid object
                    self.layout[prop] = {}

            # ### Dispatch change callbacks ###
            self._dispatch_layout_change_callbacks(delta_transform)

            # ### Layout edits no longer in process ###
            self._layout_edit_in_process = False

            # ### Call any waiting layout edit callbacks ###
            if not self._trace_edit_in_process:
                while self._waiting_edit_callbacks:
                    self._waiting_edit_callbacks.pop()()

        self._js2py_layoutDelta = None

    @observe("_js2py_restyle")
    def _handler_js2py_restyle(self, change):
        """
        Process Plotly.restyle message from the frontend
        """

        # Receive message
        # ---------------
        restyle_msg = change["new"]

        if not restyle_msg:
            self._js2py_restyle = None
            return

        style_data = restyle_msg["style_data"]
        style_traces = restyle_msg["style_traces"]
        source_view_id = restyle_msg["source_view_id"]

        # Perform restyle
        # ---------------
        self.plotly_restyle(
            restyle_data=style_data,
            trace_indexes=style_traces,
            source_view_id=source_view_id,
        )

        self._js2py_restyle = None

    @observe("_js2py_update")
    def _handler_js2py_update(self, change):
        """
        Process Plotly.update message from the frontend
        """

        # Receive message
        # ---------------
        update_msg = change["new"]

        if not update_msg:
            self._js2py_update = None
            return

        style = update_msg["style_data"]
        trace_indexes = update_msg["style_traces"]
        layout = update_msg["layout_data"]
        source_view_id = update_msg["source_view_id"]

        # Perform update
        # --------------
        self.plotly_update(
            restyle_data=style,
            relayout_data=layout,
            trace_indexes=trace_indexes,
            source_view_id=source_view_id,
        )

        self._js2py_update = None

    @observe("_js2py_relayout")
    def _handler_js2py_relayout(self, change):
        """
        Process Plotly.relayout message from the frontend
        """

        # Receive message
        # ---------------
        relayout_msg = change["new"]

        if not relayout_msg:
            self._js2py_relayout = None
            return

        relayout_data = relayout_msg["relayout_data"]
        source_view_id = relayout_msg["source_view_id"]

        if "lastInputTime" in relayout_data:
            # Remove 'lastInputTime'. Seems to be an internal plotly
            # property that is introduced for some plot types, but it is not
            # actually a property in the schema
            relayout_data.pop("lastInputTime")

        # Perform relayout
        # ----------------
        self.plotly_relayout(relayout_data=relayout_data, source_view_id=source_view_id)

        self._js2py_relayout = None

    @observe("_js2py_pointsCallback")
    def _handler_js2py_pointsCallback(self, change):
        """
        Process points callback message from the frontend
        """

        # Receive message
        # ---------------
        callback_data = change["new"]

        if not callback_data:
            self._js2py_pointsCallback = None
            return

        # Get event type
        # --------------
        event_type = callback_data["event_type"]

        # Build Selector Object
        # ---------------------
        if callback_data.get("selector", None):
            selector_data = callback_data["selector"]
            selector_type = selector_data["type"]
            selector_state = selector_data["selector_state"]
            if selector_type == "box":
                selector = BoxSelector(**selector_state)
            elif selector_type == "lasso":
                selector = LassoSelector(**selector_state)
            else:
                raise ValueError("Unsupported selector type: %s" % selector_type)
        else:
            selector = None

        # Build Input Device State Object
        # -------------------------------
        if callback_data.get("device_state", None):
            device_state_data = callback_data["device_state"]
            state = InputDeviceState(**device_state_data)
        else:
            state = None

        # Build Trace Points Dictionary
        # -----------------------------
        points_data = callback_data["points"]
        trace_points = {
            trace_ind: {
                "point_inds": [],
                "xs": [],
                "ys": [],
                "trace_name": self._data_objs[trace_ind].name,
                "trace_index": trace_ind,
            }
            for trace_ind in range(len(self._data_objs))
        }

        for x, y, point_ind, trace_ind in zip(
            points_data["xs"],
            points_data["ys"],
            points_data["point_indexes"],
            points_data["trace_indexes"],
        ):
            trace_dict = trace_points[trace_ind]
            trace_dict["xs"].append(x)
            trace_dict["ys"].append(y)
            trace_dict["point_inds"].append(point_ind)

        # Dispatch callbacks
        # ------------------
        for trace_ind, trace_points_data in trace_points.items():
            points = Points(**trace_points_data)
            trace = self.data[trace_ind]

            if event_type == "plotly_click":
                trace._dispatch_on_click(points, state)
            elif event_type == "plotly_hover":
                trace._dispatch_on_hover(points, state)
            elif event_type == "plotly_unhover":
                trace._dispatch_on_unhover(points, state)
            elif event_type == "plotly_selected":
                trace._dispatch_on_selection(points, selector)
            elif event_type == "plotly_deselect":
                trace._dispatch_on_deselect(points)

        self._js2py_pointsCallback = None

    # Display
    # -------
    def _repr_html_(self):
        """
        Customize html representation
        """
        raise NotImplementedError  # Prefer _repr_mimebundle_

    def _repr_mimebundle_(self, include=None, exclude=None, validate=True, **kwargs):
        """
        Return mimebundle corresponding to default renderer.
        """
        display_jupyter_version_warnings()

        # Widget layout and data need to be set here in case there are
        # changes made to the figure after the widget is created but before
        # the cell is run.
        self._widget_layout = deepcopy(self._layout_obj._props)
        self._widget_data = deepcopy(self._data)
        return {
            "application/vnd.jupyter.widget-view+json": {
                "version_major": 2,
                "version_minor": 0,
                "model_id": self._model_id,
            },
        }

    def _ipython_display_(self):
        """
        Handle rich display of figures in ipython contexts
        """
        raise NotImplementedError  # Prefer _repr_mimebundle_

    # Callbacks
    # ---------
    def on_edits_completed(self, fn):
        """
        Register a function to be called after all pending trace and layout
        edit operations have completed

        If there are no pending edit operations then function is called
        immediately

        Parameters
        ----------
        fn : callable
            Function of zero arguments to be called when all pending edit
            operations have completed
        """
        if self._layout_edit_in_process or self._trace_edit_in_process:
            self._waiting_edit_callbacks.append(fn)
        else:
            fn()

    # Validate No Frames
    # ------------------
    @property
    def frames(self):
        # Note: This property getter is identical to that of the superclass,
        # but it must be included here because we're overriding the setter
        # below.
        return self._frame_objs

    @frames.setter
    def frames(self, new_frames):
        if new_frames:
            BaseFigureWidget._display_frames_error()

    @staticmethod
    def _display_frames_error():
        """
        Display an informative error when user attempts to set frames on a
        FigureWidget

        Raises
        ------
        ValueError
            always
        """
        msg = """
Frames are not supported by the plotly.graph_objs.FigureWidget class.
Note: Frames are supported by the plotly.graph_objs.Figure class"""
        raise ValueError(msg)

    # Static Helpers
    # --------------
    @staticmethod
    def _remove_overlapping_props(input_data, delta_data, prop_path=()):
        """
        Remove properties in input_data that are also in delta_data, and do so
        recursively.

        Exception: Never remove 'uid' from input_data, this property is used
        to align traces

        Parameters
        ----------
        input_data : dict|list
        delta_data : dict|list

        Returns
        -------
        list[tuple[str|int]]
            List of removed property path tuples
        """

        # Initialize removed
        # ------------------
        # This is the list of path tuples to the properties that were
        # removed from input_data
        removed = []

        # Handle dict
        # -----------
        if isinstance(input_data, dict):
            assert isinstance(delta_data, dict)

            for p, delta_val in delta_data.items():
                if isinstance(delta_val, dict) or BaseFigure._is_dict_list(delta_val):
                    if p in input_data:
                        # ### Recurse ###
                        input_val = input_data[p]
                        recur_prop_path = prop_path + (p,)
                        recur_removed = BaseFigureWidget._remove_overlapping_props(
                            input_val, delta_val, recur_prop_path
                        )
                        removed.extend(recur_removed)

                        # Check whether the last property in input_val
                        # has been removed. If so, remove it entirely
                        if not input_val:
                            input_data.pop(p)
                            removed.append(recur_prop_path)

                elif p in input_data and p != "uid":
                    # ### Remove property ###
                    input_data.pop(p)
                    removed.append(prop_path + (p,))

        # Handle list
        # -----------
        elif isinstance(input_data, list):
            assert isinstance(delta_data, list)

            for i, delta_val in enumerate(delta_data):
                if i >= len(input_data):
                    break

                input_val = input_data[i]
                if (
                    input_val is not None
                    and isinstance(delta_val, dict)
                    or BaseFigure._is_dict_list(delta_val)
                ):
                    # ### Recurse ###
                    recur_prop_path = prop_path + (i,)
                    recur_removed = BaseFigureWidget._remove_overlapping_props(
                        input_val, delta_val, recur_prop_path
                    )

                    removed.extend(recur_removed)

        return removed

    @staticmethod
    def _transform_data(to_data, from_data, should_remove=True, relayout_path=()):
        """
        Transform to_data into from_data and return relayout-style
        description of the transformation

        Parameters
        ----------
        to_data : dict|list
        from_data : dict|list

        Returns
        -------
        dict
            relayout-style description of the transformation
        """

        # Initialize relayout data
        # ------------------------
        relayout_data = {}

        # Handle dict
        # -----------
        if isinstance(to_data, dict):
            # ### Validate from_data ###
            if not isinstance(from_data, dict):
                raise ValueError(
                    "Mismatched data types: {to_dict} {from_data}".format(
                        to_dict=to_data, from_data=from_data
                    )
                )

            # ### Add/modify properties ###
            # Loop over props/vals
            for from_prop, from_val in from_data.items():
                # #### Handle compound vals recursively ####
                if isinstance(from_val, dict) or BaseFigure._is_dict_list(from_val):
                    # ##### Init property value if needed #####
                    if from_prop not in to_data:
                        to_data[from_prop] = {} if isinstance(from_val, dict) else []

                    # ##### Transform property val recursively #####
                    input_val = to_data[from_prop]
                    relayout_data.update(
                        BaseFigureWidget._transform_data(
                            input_val,
                            from_val,
                            should_remove=should_remove,
                            relayout_path=relayout_path + (from_prop,),
                        )
                    )

                # #### Handle simple vals directly ####
                else:
                    if from_prop not in to_data or not BasePlotlyType._vals_equal(
                        to_data[from_prop], from_val
                    ):
                        to_data[from_prop] = from_val
                        relayout_path_prop = relayout_path + (from_prop,)
                        relayout_data[relayout_path_prop] = from_val

            # ### Remove properties ###
            if should_remove:
                for remove_prop in set(to_data.keys()).difference(
                    set(from_data.keys())
                ):
                    to_data.pop(remove_prop)

        # Handle list
        # -----------
        elif isinstance(to_data, list):
            # ### Validate from_data ###
            if not isinstance(from_data, list):
                raise ValueError(
                    "Mismatched data types: to_data: {to_data} {from_data}".format(
                        to_data=to_data, from_data=from_data
                    )
                )

            # ### Add/modify properties ###
            # Loop over indexes / elements
            for i, from_val in enumerate(from_data):
                # #### Initialize element if needed ####
                if i >= len(to_data):
                    to_data.append(None)
                input_val = to_data[i]

                # #### Handle compound element recursively ####
                if input_val is not None and (
                    isinstance(from_val, dict) or BaseFigure._is_dict_list(from_val)
                ):
                    relayout_data.update(
                        BaseFigureWidget._transform_data(
                            input_val,
                            from_val,
                            should_remove=should_remove,
                            relayout_path=relayout_path + (i,),
                        )
                    )

                # #### Handle simple elements directly ####
                else:
                    if not BasePlotlyType._vals_equal(to_data[i], from_val):
                        to_data[i] = from_val
                        relayout_data[relayout_path + (i,)] = from_val

        return relayout_data
