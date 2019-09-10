.. raw:: html

   <script>
   window.addEventListener('DOMContentLoaded', function() {
        (function($) {
   //Function to make the index toctree collapsible
   $(function () {
       $('div.body .toctree-l2')
           .click(function(event){
               if (event.target.tagName.toLowerCase() != "a") {
                   if ($(this).children('ul').length > 0) {
                        $(this).attr('data-content',
                            (!$(this).children('ul').is(':hidden')) ? '\u25ba' : '\u25bc');
                       $(this).children('ul').toggle();
                   }
                   return true; //Makes links clickable
               }
           })
           .mousedown(function(event){ return false; }) //Firefox highlighting fix
           .children('ul').hide();
       // Initialize the values
       $('div.body li.toctree-l2:not(:has(ul))').attr('data-content', '-');
       $('div.body li.toctree-l2:has(ul)').attr('data-content', '\u25ba');
       $('div.body li.toctree-l2:has(ul)').css('cursor', 'pointer');

       $('div.body .toctree-l2').hover(
           function () {
               if ($(this).children('ul').length > 0) {
                   $(this).css('background-color', '#e5e5e5').children('ul').css('background-color', '#F0F0F0');
                   $(this).attr('data-content',
                       (!$(this).children('ul').is(':hidden')) ? '\u25bc' : '\u25ba');
               }
               else {
                   $(this).css('background-color', '#F9F9F9');
               }
           },
           function () {
               $(this).css('background-color', 'white').children('ul').css('background-color', 'white');
               if ($(this).children('ul').length > 0) {
                   $(this).attr('data-content',
                       (!$(this).children('ul').is(':hidden')) ? '\u25bc' : '\u25ba');
               }
           }
       );
   });
        })(jQuery);
    });
   </script>

  <style type="text/css">
    div.body li, div.body ul {
        transition-duration: 0.2s;
    }

    div.body li.toctree-l1 {
        padding: 5px 0 0;
        list-style-type: none;
        font-size: 150%;
        background-color: #f2f2f2;
        font-weight: normal;
        color: #20435c;
        margin-left: 0;
        margin-bottom: 1.2em;
        font-weight: bold;
        }

    div.body li.toctree-l1 a {
        color: #314F64;
    }

    div.body li.toctree-l1 > a {
        margin-left: 0.75rem;
    }

    div.body li.toctree-l2 {
        padding: 0.25em 0 0.25em 0 ;
        list-style-type: none;
        background-color: #FFFFFF;
        font-size: 85% ;
        font-weight: normal;
        margin-left: 0;
    }

    div.body li.toctree-l2 ul {
        padding-left: 40px ;
    }

    div.body li.toctree-l2:before {
        content: attr(data-content);
        font-size: 1rem;
        color: #777;
        display: inline-block;
        width: 1.5rem;
    }

    div.body li.toctree-l3 {
        font-size: 88% ;
        list-style-type: square;
        font-weight: normal;
        margin-left: 0;
    }

    div.body li.toctree-l4 {
        font-size: 93% ;
        list-style-type: circle;
        font-weight: normal;
        margin-left: 0;
    }

    div.body div.topic li.toctree-l1 {
        font-size: 100% ;
        font-weight: bold;
        background-color: transparent;
        margin-bottom: 0;
        margin-left: 1.5em;
        display:inline;
    }

    div.body div.topic p {
        font-size: 90% ;
        margin: 0.4ex;
    }

    div.body div.topic p.topic-title {
        display:inline;
        font-size: 100% ;
        margin-bottom: 0;
    }
  </style>


