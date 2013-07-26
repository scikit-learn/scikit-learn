.. title:: User guide: contents

..
    We are putting the title as a raw HTML so that it doesn't appear in
    the contents

.. raw:: html

    <h1>User guide: contents</h1>

    <script>
    //Function to make the index toctree collapsible
        $(function () {
                $('.toctree-l2')
                    .click(function(event){
                        if (event.target.tagName.toLowerCase() != "a") {
                    if ($(this).children('ul').length > 0) {
                                $(this).css('list-style-image',
                                (!$(this).children('ul').is(':hidden')) ? 'url(_static/plusBoxHighlight.png)' : 'url(_static/minBoxHighlight.png)');
                                $(this).children('ul').toggle("slow");
                            }
                            return true; //Makes links clickable
                        }
            })
            .mousedown(function(event){ return false; }) //Firefox highlighting fix
                    .css({cursor:'pointer', 'list-style-image':'url(_static/plusBox.png)'})
                    .children('ul').hide();
                $('ul li ul li:not(:has(ul))').css({cursor:'default', 'list-style-image':'url(_static/noneBox.png)'});
                $('ul li ul').css('margin-left', '0px');
            $('.toctree-l3').css({cursor:'default', 'list-style-image':'url(_static/noneBox.png)'});
            $('.toctree-l2').hover(
                function () {
                if ($(this).children('ul').length > 0) {
                    $(this).css('background-color', '#D0D0D0').children('ul').css('background-color', '#F0F0F0');
                    $(this).css('list-style-image',
                                (!$(this).children('ul').is(':hidden')) ? 'url(_static/minBoxHighlight.png)' : 'url(_static/plusBoxHighlight.png)');
                }
                else {
                    $(this).css('background-color', '#F9F9F9');
                }
                    },
                    function () {
                        $(this).css('background-color', 'white').children('ul').css('background-color', 'white');
                if ($(this).children('ul').length > 0) {
                    $(this).css('list-style-image',
                                (!$(this).children('ul').is(':hidden')) ? 'url(_static/minBox.png)' : 'url(_static/plusBox.png)');
                }
                    }
                );
        });
    </script>

.. _user_guide:

.. include:: includes/big_toc_css.rst

.. toctree::
   :numbered:

   supervised_learning.rst
   unsupervised_learning.rst
   model_selection.rst
   data_transforms.rst
   Dataset loading utilities <datasets/index.rst>
