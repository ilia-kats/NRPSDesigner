jQuery(document).ready(function() {
			
	if (jQuery('#user_tab_profile').length == 1) {
		jQuery('#user_tab_profile').button({
			icons:{primary:'ui-icon-person',
				   secondary:'ui-icon-triangle-1-s'}
		});
		jQuery('#user_tab_inbox').button({
			icons:{primary:'ui-icon-mail-closed'}
		});
	} else {
		jQuery('#user_tab_login').button({
			icons:{primary:'ui-icon-tag'}
		});
		jQuery('#user_tab_register').button({
			icons:{primary:'ui-icon-key'}
		});
	}
	jQuery('#user_tab_home').button({
		icons:{primary:'ui-icon-home',
			   secondary:'ui-icon-triangle-1-s'}
	});
	jQuery('#user_tab_help').button({
		icons:{primary:'ui-icon-help',
			  }
	});
	jQuery('#user_tab_tools').button({
		icons:{primary:'ui-icon-wrench',
			   secondary:'ui-icon-triangle-1-s'}
	});
	
	jQuery('.user_tab').addClass('ui-corner-top').removeClass('ui-corner-all');

//	$('span.drop').mouseleave(function () {$('ul', this).slideUp("fast") });
	
//	$('span.drop a.user_tab').mouseenter(function () {$('~ ul', this).slideDown('fast')});

	jQuery('span.drop').mouseleave(function () {jQuery('ul', this).hide()});
	jQuery('span.drop a.user_tab').mouseenter(function() {jQuery('~ ul', this).show()});
	
	jQuery('span.drop ul li a').button();
	
	jQuery('input[id^="id_captcha"]').addClass('captcha');
});
  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-18410656-4']);
  _gaq.push(['_trackPageview']);
