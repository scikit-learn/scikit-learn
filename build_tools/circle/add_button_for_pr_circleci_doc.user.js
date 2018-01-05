// ==UserScript==
// @name        see PR doc on CircleCI
// @namespace   scikit-learn.org
// @include     /https://github.com/scikit-learn/scikit-learn/pull/[0-9]+[^/]*$/
// @grant       none
// ==/UserScript==

(function(){
    'use strict';
    // console.log('Running user script');

    window.addEventListener('load', () => {
        addButton('See CircleCI doc for this PR', gotoCirclePage);
    });

    function addButton(text, onclick, cssObj) {
        let headerActionElement = document.getElementsByClassName('gh-header-actions')[0];
        // console.log('headerActionElement', headerActionElement);

        cssObj = cssObj || {};
        let button = document.createElement('button'), btnStyle = button.style;
        headerActionElement.appendChild(button);
        button.innerHTML = text;
        button.onclick = onclick;
        button.classList = "btn btn-sm";
        button.type = "button";
        Object.keys(cssObj).forEach(key => btnStyle[key] = cssObj[key]);
        return button;
    }

    function findElementsBySelectorAndText(selector, text) {
        var elements = document.querySelectorAll(selector);
        return Array.prototype.filter.call(elements, function(element){
            return RegExp(text).test(element.textContent);
        });
    }

    function gotoCirclePage() {
        var circleElement;
        var useCircleWorkflow;
        circleElement = findElementsBySelectorAndText('.branch-action .merge-status-item',
                                                      'ci/circle.+python3')[0];

        if (circleElement) {
            useCircleWorkflow = true;
        } else {
            circleElement = findElementsBySelectorAndText('.branch-action .merge-status-item',
                                                          'ci/circle')[0];
            useCircleWorkflow = false;
        }
        // console.log('circleElement', circleElement);

        var match = /circleci.com\/.*?([0-9]+)\?/.exec(circleElement.innerHTML);
        // console.log('match', match);

        if (match) {
            var circleBuildNumber = match[1];
            var docURLPart;
            if (useCircleWorkflow) {
                docURLPart = '-843222-gh.circle-artifacts.com/0/doc';
            } else {
                docURLPart = '-843222-gh.circle-artifacts.com/0/home/ubuntu/scikit-learn/doc/_build/html/stable';
            }
            var docURL = 'https://' + circleBuildNumber + docURLPart + '/_changed.html';
            // console.log('docURL', docURL);
            window.open(docURL);
        }
    }
}());
