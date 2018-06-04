var comm_websocket_adapter = function(comm) {
    // Create a "websocket"-like object which calls the given IPython comm
    // object with the appropriate methods. Currently this is a non binary
    // socket, so there is still some room for performance tuning.
    var ws = {};

    ws.close = function() {
        comm.close()
    };
    ws.send = function(m) {
        //console.log('sending', m);
        comm.send(m);
    };
    // Register the callback with on_msg.
    comm.on_msg(function(msg) {
        //console.log('receiving', msg['content']['data'], msg);
        // Pass the mpl event to the overridden (by mpl) onmessage function.
        ws.onmessage(msg['content']['data'])
    });
    return ws;
}

mpl.mpl_figure_comm = function(comm, msg) {
    // This is the function which gets called when the mpl process
    // starts-up an IPython Comm through the "matplotlib" channel.

    var id = msg.content.data.id;
    // Get hold of the div created by the display call when the Comm
    // socket was opened in Python.
    var element = $("#" + id);
    var ws_proxy = comm_websocket_adapter(comm)

    function ondownload(figure, format) {
        window.open(figure.imageObj.src);
    }

    var fig = new mpl.figure(id, ws_proxy,
                           ondownload,
                           element.get(0));

    // Call onopen now - mpl needs it, as it is assuming we've passed it a real
    // web socket which is closed, not our websocket->open comm proxy.
    ws_proxy.onopen();

    fig.parent_element = element.get(0);
    fig.cell_info = mpl.find_output_cell("<div id='" + id + "'></div>");
    if (!fig.cell_info) {
        console.error("Failed to find cell for figure", id, fig);
        return;
    }

    var output_index = fig.cell_info[2]
    var cell = fig.cell_info[0];

};

mpl.figure.prototype.handle_close = function(fig, msg) {
    var width = fig.canvas.width/mpl.ratio
    fig.root.unbind('remove')

    // Update the output cell to use the data from the current canvas.
    fig.push_to_output();
    var dataURL = fig.canvas.toDataURL();
    // Re-enable the keyboard manager in IPython - without this line, in FF,
    // the notebook keyboard shortcuts fail.
    IPython.keyboard_manager.enable()
    $(fig.parent_element).html('<img src="' + dataURL + '" width="' + width + '">');
    fig.close_ws(fig, msg);
}

mpl.figure.prototype.close_ws = function(fig, msg){
    fig.send_message('closing', msg);
    // fig.ws.close()
}

mpl.figure.prototype.push_to_output = function(remove_interactive) {
    // Turn the data on the canvas into data in the output cell.
    var width = this.canvas.width/mpl.ratio
    var dataURL = this.canvas.toDataURL();
    this.cell_info[1]['text/html'] = '<img src="' + dataURL + '" width="' + width + '">';
}

mpl.figure.prototype.updated_canvas_event = function() {
    // Tell IPython that the notebook contents must change.
    IPython.notebook.set_dirty(true);
    this.send_message("ack", {});
    var fig = this;
    // Wait a second, then push the new image to the DOM so
    // that it is saved nicely (might be nice to debounce this).
    setTimeout(function () { fig.push_to_output() }, 1000);
}

mpl.figure.prototype._init_toolbar = function() {
    var fig = this;

    var nav_element = $('<div/>')
    nav_element.attr('style', 'width: 100%');
    this.root.append(nav_element);

    // Define a callback function for later on.
    function toolbar_event(event) {
        return fig.toolbar_button_onclick(event['data']);
    }
    function toolbar_mouse_event(event) {
        return fig.toolbar_button_onmouseover(event['data']);
    }

    for(var toolbar_ind in mpl.toolbar_items){
        var name = mpl.toolbar_items[toolbar_ind][0];
        var tooltip = mpl.toolbar_items[toolbar_ind][1];
        var image = mpl.toolbar_items[toolbar_ind][2];
        var method_name = mpl.toolbar_items[toolbar_ind][3];

        if (!name) { continue; };

        var button = $('<button class="btn btn-default" href="#" title="' + name + '"><i class="fa ' + image + ' fa-lg"></i></button>');
        button.click(method_name, toolbar_event);
        button.mouseover(tooltip, toolbar_mouse_event);
        nav_element.append(button);
    }

    // Add the status bar.
    var status_bar = $('<span class="mpl-message" style="text-align:right; float: right;"/>');
    nav_element.append(status_bar);
    this.message = status_bar[0];

    // Add the close button to the window.
    var buttongrp = $('<div class="btn-group inline pull-right"></div>');
    var button = $('<button class="btn btn-mini btn-primary" href="#" title="Stop Interaction"><i class="fa fa-power-off icon-remove icon-large"></i></button>');
    button.click(function (evt) { fig.handle_close(fig, {}); } );
    button.mouseover('Stop Interaction', toolbar_mouse_event);
    buttongrp.append(button);
    var titlebar = this.root.find($('.ui-dialog-titlebar'));
    titlebar.prepend(buttongrp);
}

mpl.figure.prototype._root_extra_style = function(el){
    var fig = this
    el.on("remove", function(){
	fig.close_ws(fig, {});
    });
}

mpl.figure.prototype._canvas_extra_style = function(el){
    // this is important to make the div 'focusable
    el.attr('tabindex', 0)
    // reach out to IPython and tell the keyboard manager to turn it's self
    // off when our div gets focus

    // location in version 3
    if (IPython.notebook.keyboard_manager) {
        IPython.notebook.keyboard_manager.register_events(el);
    }
    else {
        // location in version 2
        IPython.keyboard_manager.register_events(el);
    }

}

mpl.figure.prototype._key_event_extra = function(event, name) {
    var manager = IPython.notebook.keyboard_manager;
    if (!manager)
        manager = IPython.keyboard_manager;

    // Check for shift+enter
    if (event.shiftKey && event.which == 13) {
        this.canvas_div.blur();
        event.shiftKey = false;
        // Send a "J" for go to next cell
        event.which = 74;
        event.keyCode = 74;
        manager.command_mode();
        manager.handle_keydown(event);
    }
}

mpl.figure.prototype.handle_save = function(fig, msg) {
    fig.ondownload(fig, null);
}


mpl.find_output_cell = function(html_output) {
    // Return the cell and output element which can be found *uniquely* in the notebook.
    // Note - this is a bit hacky, but it is done because the "notebook_saving.Notebook"
    // IPython event is triggered only after the cells have been serialised, which for
    // our purposes (turning an active figure into a static one), is too late.
    var cells = IPython.notebook.get_cells();
    var ncells = cells.length;
    for (var i=0; i<ncells; i++) {
        var cell = cells[i];
        if (cell.cell_type === 'code'){
            for (var j=0; j<cell.output_area.outputs.length; j++) {
                var data = cell.output_area.outputs[j];
                if (data.data) {
                    // IPython >= 3 moved mimebundle to data attribute of output
                    data = data.data;
                }
                if (data['text/html'] == html_output) {
                    return [cell, data, j];
                }
            }
        }
    }
}

// Register the function which deals with the matplotlib target/channel.
// The kernel may be null if the page has been refreshed.
if (IPython.notebook.kernel != null) {
    IPython.notebook.kernel.comm_manager.register_target('matplotlib', mpl.mpl_figure_comm);
}
