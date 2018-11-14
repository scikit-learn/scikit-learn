/*
MIT License

Copyright (c) 2018 Chris Holdgraf

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

// Localization support
const messages = {
    'en': {
      'copy': 'Copy   except prompts and output',
      'copy_to_clipboard': 'Copy to clipboard',
      'copy_success': 'Copied except prompts and output',
      'copy_failure': 'Failed to copy',
    }
  }

  let locale = 'en'
  if( document.documentElement.lang !== undefined
      && messages[document.documentElement.lang] !== undefined ) {
    locale = document.documentElement.lang
  }

  /**
   * Set up copy/paste for code blocks
   */

  const runWhenDOMLoaded = cb => {
    if (document.readyState != 'loading') {
      cb()
    } else if (document.addEventListener) {
      document.addEventListener('DOMContentLoaded', cb)
    } else {
      document.attachEvent('onreadystatechange', function() {
        if (document.readyState == 'complete') cb()
      })
    }
  }

  const codeCellId = index => `codecell${index}`

  // Clears selected text since ClipboardJS will select the text when copying
  const clearSelection = () => {
    if (window.getSelection) {
      window.getSelection().removeAllRanges()
    } else if (document.selection) {
      document.selection.empty()
    }
  }

  // Changes tooltip text for two seconds, then changes it back
  const temporarilyChangeTooltip = (el, newText) => {
    const oldText = el.getAttribute('data-tooltip')
    el.setAttribute('data-tooltip', newText)
    setTimeout(() => el.setAttribute('data-tooltip', oldText), 2000)
  }

  const addCopyButtonToCodeCells = () => {
    // If ClipboardJS hasn't loaded, wait a bit and try again. This
    // happens because we load ClipboardJS asynchronously.
    if (window.ClipboardJS === undefined) {
      setTimeout(addCopyButtonToCodeCells, 250)
      return
    }

  const codeCells = document.querySelectorAll('div.highlight pre')
    codeCells.forEach((codeCell, index) => {
      const id = codeCellId(index)
      codeCell.setAttribute('id', id)

      const clipboardButton = id =>
      `<a class="copybtn o-tooltip--right" data-tooltip="${messages[locale]['copy']}" data-clipboard-target="#${id}">
        <img src="https://gitcdn.xyz/repo/choldgraf/sphinx-copybutton/master/sphinx_copybutton/_static/copy-button.svg" alt="${messages[locale]['copy_to_clipboard']}">
      </a>`
      codeCell.insertAdjacentHTML('afterend', clipboardButton(id))
    })

  var class_to_remove = ['gp','go','gt'];
  function toggle_prompt_output(codeblock,class_to_remove,state){
    var subs = codeblock.getElementsByClassName(class_to_remove);
    for(var i = 0; i < subs.length; i++){
      var a = subs[i];
      a.style.display = state;
      }
    }

  const clipboard = new ClipboardJS('.copybtn', {
    text: function(trigger) {
    const query = trigger.getAttribute('data-clipboard-target');
    var htmlBlock = document.querySelectorAll('div.highlight pre'+query)[0]
    //Hide Outputs and Prompts
    for (i = 0; i < class_to_remove.length; i++) {
      toggle_prompt_output(htmlBlock, class_to_remove[i], 'none')
    }
    return htmlBlock.innerText
      }
    });

  clipboard.on('success', event => {
    const query = event.trigger.getAttribute('data-clipboard-target');
    var htmlBlock = document.querySelectorAll('div.highlight pre'+query)[0]
    //Show Outputs and Prompts
    for (i = 0; i < class_to_remove.length; i++) {
      toggle_prompt_output(htmlBlock, class_to_remove[i], 'inline')
    }

    clearSelection()
    temporarilyChangeTooltip(event.trigger, messages[locale]['copy_success'])
  })

  clipboard.on('error', event => {
    temporarilyChangeTooltip(event.trigger, messages[locale]['copy_failure'])
  })
  }

  runWhenDOMLoaded(addCopyButtonToCodeCells)
