// Copyright 2013 PSF. Licensed under the PYTHON SOFTWARE FOUNDATION LICENSE VERSION 2
// File originates from the cpython source found in Doc/tools/sphinxext/static/version_switch.js

(function() {
  'use strict';

  var all_versions = {
    'latest': '0.19',
    'v0.18': '0.18',
    'v0.17': '0.17',
    'v0.16': '0.16',
    'v0.15': '0.15',
    'v0.14': '0.14',
    'v0.13': '0.13',
    'v0.12': '0.12',
    'v0.11': '0.11',
    'v0.10': '0.10',
    'v0.9': '0.9',
    'v0.8': '0.8',
    'v0.7': '0.7',
    'v0.6': '0.6',
    'v0.5': '0.5',
    'v0.4': '0.4'
  };

  function build_select(current_version, current_release) {
    var buf = ['<select>'];

    $.each(all_versions, function(path, version) {
      buf.push('<option value="' + path + '"');
      if (version == current_version)
        buf.push(' selected="selected">' + current_release + '</option>');
      else
        buf.push('>' + version + '</option>');
    });

    buf.push('</select>');
    return buf.join('');
  }

  function patch_url(url, new_version) {
    var url_re = /scitools\.org\.uk\/cartopy\/docs\/(latest|(v\d+\.\d+))\//,
        new_url = url.replace(url_re, 'scitools.org.uk/cartopy/docs/' + new_version + '/');

    if (new_url == url && !new_url.match(url_re)) {
      alert("The version switcher only functions on the scitools.org.uk domain.");
    }
    return new_url;
  }

  function on_switch() {
    var selected = $(this).children('option:selected').attr('value');

    var url = window.location.href,
        new_url = patch_url(url, selected);

    if (new_url != url) {
      // check beforehand if url exists, else redirect to version's start page
      $.ajax({
        url: new_url,
        success: function() {
           window.location.href = new_url;
        },
        error: function() {
           window.location.href = 'https://scitools.org.uk/cartopy/docs/' + selected;
        }
      });
    }
  }

  $(document).ready(function() {
    var release = DOCUMENTATION_OPTIONS.VERSION;
    // Take the first 2 parts of the release (e.g. "0.6.0" -> "0.6")
    var version = release.split('.').slice(0, 2).join('.');
    var select = build_select(version, release);

    var index_li = $('li.right:contains("index")');
    index_li.append('|&nbsp;'); 
    index_li.before('<li class="version_switcher right">' + select + '</li>');
    $('.version_switcher select').bind('change', on_switch);
  });
})();
