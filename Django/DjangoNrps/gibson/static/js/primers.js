function makeHttpObject() {
  try {return new XMLHttpRequest();}
  catch (error) {}
  try {return new ActiveXObject("Msxml2.XMLHTTP");}
  catch (error) {}
  try {return new ActiveXObject("Microsoft.XMLHTTP");}
  catch (error) {}

  throw new Error("Could not create HTTP request object.");
}

var validate = function(input) {
	re = /^\d+(?:\.\d*)?jQuery/;
	if(re.exec(input.value) === null){
		jQuery(input.parentNode).addClass("errorcell");
		return false;
	}else{
		jQuery(input.parentNode).removeClass("errorcell");
		return true;
	}
}

var Mix = function () {
	this.enz_d = jQuery('#enzyme_desired')[0];
	this.enz_s = jQuery('#enzyme_stock')[0];
	this.buff_d = jQuery('#buffer_desired')[0];
	this.buff_s = jQuery('#buffer_stock')[0];
	this.dntp_d = jQuery('#dntp_desired')[0];
	this.dntp_s = jQuery('#dntp_stock')[0];
	this.prim_d = jQuery('#primer_desired')[0];
	this.temp_d = jQuery('#template_desired')[0];
	this.repeats = jQuery('#repeats')[0];
	this.error = jQuery('#error_margin')[0];
	this.vol_e = jQuery('#volume_each')[0];
	var mixform = jQuery('form.mixform');
	this.errors = 0;
	this.warnings = 0;
	this.fragments = []
	for (f in mixform){
		if (mixform[f].constructor == HTMLFormElement){
			this.fragments.push(new MixFragment(this, mixform[f]));
		}
	}
	this.go();
}

Mix.prototype.validate = function () {
	this.errors = 0;
	for (x in this){
		if (this[x].constructor == HTMLInputElement){
			this.errors += !validate(this[x]);
		}
	}
}

Mix.prototype.m = function () {
	return this.repeats.value * (1 + this.error.value/100);
}

Mix.prototype.go = function () {
	this.warnings = 0;
	this.validate();
	for (f in this.fragments) {
		this.fragments[f].validate();
	}
	if (this.errors > 0) {
		return false;
	}
	for (f in this.fragments) {
		this.fragments[f].go();
		this.fragments[f].warn();
	}
	if (this.warnings > 0) {
		jQuery('#warning').show();
	} else {
		jQuery('#warning').hide();
	}
}


var MixFragment = function (mix, form) {
	this.mix = mix;
	this.prim_t_s = form.elements['primer_top_stock'];
	this.prim_b_s = form.elements['primer_bottom_stock'];
	this.temp_s = form.elements['template_stock'];
	this.enz_v = form.elements['enzyme_vol'];
	this.buff_v = form.elements['buffer_vol'];
	this.dntp_v = form.elements['dntp_vol'];
	this.prim_t_v = form.elements['primer_top_vol'];
	this.prim_b_v = form.elements['primer_bottom_vol'];
	this.temp_v = form.elements['template_vol'];
	this.water_v = form.elements['water_vol'];
	this.total_v = form.elements['total_vol'];
}

MixFragment.prototype.go = function () {
	this.enz_v.value = (this.mix.m()*this.mix.enz_d.value/this.mix.enz_s.value).toFixed(1);
	this.buff_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.buff_d.value/this.mix.buff_s.value).toFixed(1);
	this.dntp_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.dntp_d.value/this.mix.dntp_s.value).toFixed(1);
	this.temp_v.value = (this.mix.m()*this.mix.temp_d.value/this.temp_s.value).toFixed(1);
	this.prim_t_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.prim_d.value/this.prim_t_s.value).toFixed(1);
	this.prim_b_v.value = (this.mix.vol_e.value*this.mix.m()*this.mix.prim_d.value/this.prim_b_s.value).toFixed(1);
	this.water_v.value = (this.mix.m()*this.mix.vol_e.value - (parseFloat(this.enz_v.value) + parseFloat(this.buff_v.value) + parseFloat(this.dntp_v.value) + parseFloat(this.temp_v.value) + parseFloat(this.prim_t_v.value) + parseFloat(this.prim_b_v.value))).toFixed(1);
	this.total_v.value = (this.mix.m()*this.mix.vol_e.value).toFixed(1);
}

MixFragment.prototype.validate = function () {
	for (x in this){
		if (this[x].constructor == HTMLInputElement && !this[x].readOnly){
			this.mix.errors += !validate(this[x]);
		}
	}
}

MixFragment.prototype.warn = function () {
	for (x in this){
		if (this[x].constructor == HTMLInputElement && this[x].readOnly){
			if (parseFloat(this[x].value) < 1 ){
				jQuery(this[x].parentNode).addClass("warningcell");
				this.mix.warnings += 1;
			} else {
				jQuery(this[x].parentNode).removeClass("warningcell");
			}
		}
	}
}

jQuery('document').ready(function () {
	jQuery('#protocol-accordion').accordion({
		collapsible:true,
	});
	jQuery('#warning').hide();
	jQuery('.mix').change(function(){
		validate(this);
		if (jQuery('.errorcell').size() > 0){
			jQuery('button#download_protocol').button('disable');
		} else {
			jQuery('button#download_protocol').button('enable');
		}
	});
});
