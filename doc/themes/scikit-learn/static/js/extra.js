// Miscellaneous enhancements to doc display


$(document).ready(function() {
	/*** Add permalink buttons next to glossary terms ***/

	$('dl.glossary > dt[id]').append(function() {
		return ('<a class="headerlink" href="#' +
			    this.getAttribute('id') +
			    '" title="Permalink to this term">¶</a>');
	})
});
