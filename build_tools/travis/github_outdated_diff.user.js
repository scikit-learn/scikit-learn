// ==UserScript==
// @name        github expand outdated diff comments
// @namespace   scikit-learn.org
// @include     https://github.com/*/pull/*
// @version     1
// @grant       none
// ==/UserScript==

(function(){
    opened = false;

    document.addEventListener('keydown', function(e) {
        // console.log('Pressed: ' + e.keyCode);
        // console.log('shift: ' + e.shiftKey);
        // console.log('ctrl: ' + e.ctrlKey);
        // console.log('alt: ', e.altKey);
        // console.log('meta: ' + e.metaKey);
        if (e.keyCode == 80 && !e.shiftKey && e.ctrlKey && e.altKey && !e.metaKey) {
            // console.log('Pressed alt-ctrl-p')
            outdated_diff_elements = document.getElementsByClassName("outdated-comment");
            // supporting old github PRs (before "Start Review" feature)
            legacy_outdated_diff_elements = document.getElementsByClassName("outdated-diff-comment-container");
            // transform HTMLCollection objects into arrays and concatenate
            outdated_diff_elements = Array.concat(Array.from(outdated_diff_elements),
                                                  Array.from(legacy_outdated_diff_elements))
            console.log('outdated_diff_elements.length: ' + outdated_diff_elements.length);
            for (i = 0; i < outdated_diff_elements.length; i++) {
                element = outdated_diff_elements[i];
                if (!opened) {
                   element.classList.add('open');
               }
               else {
                   element.classList.remove('open');
               }
            }
            opened = !opened;
            // console.log('Opened: ' + opened);
        }
    }, false);
})();
