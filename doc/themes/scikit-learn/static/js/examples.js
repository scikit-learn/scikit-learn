/* Interactive component of example galleries */

function animateClone(e){
  var position;
  position = $(this).position();
  var clone = $(this).closest('.thumbnailContainer').find('.clonedItem');
  var clone_fig = clone.find('.figure');
  clone.css("left", position.left - 70).css("top", position.top - 70).css("position", "absolute").css("z-index", 1000).css("background-color", "transparent");

  var cloneImg = clone_fig.find('img');

  clone.show();
  clone.animate({
        height: "270px",
        width: "320px"
    }, 0
  );
  cloneImg.css({
        'max-height': "200px",
        'max-width': "280px"
  });
  cloneImg.animate({
        height: "200px",
        width: "280px"
    }, 0
   );
  clone_fig.css({
       'margin-top': '20px',
  });
  clone_fig.show();
  clone.find('p').css("display", "block");
  clone_fig.css({
       height: "240",
       width: "305px"
  });
  cloneP_height = clone.find('p.caption').height();
  clone_fig.animate({
       height: (200 + cloneP_height)
   }, 0
  );

  clone.bind("mouseleave", function(e){
      clone.animate({
          height: "100px",
          width: "150px"
      }, 10, function(){$(this).hide();});
      clone_fig.animate({
          height: "100px",
          width: "150px"
      }, 10, function(){$(this).hide();});
  });
} //end animateClone()


$(window).load(function () {
    $(".thumbnailContainer .figure").css("z-index", 1);

    $(".docstringWrapper").each(function(i, obj){
        var clone;
        var $obj = $(obj);
        clone = $obj.clone();
        clone.addClass("clonedItem");
        clone.appendTo($obj.closest(".thumbnailContainer"));
        clone.hide();
        $obj.bind("mouseenter", animateClone);
    }); // end each
}); // end

