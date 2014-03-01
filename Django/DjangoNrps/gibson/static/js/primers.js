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
	re = /^\d+(?:\.\d*)?/;
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
	this.prim_t_s = form.elements['p_t_c'];
	this.prim_b_s = form.elements['p_b_c'];
	this.temp_s = form.elements['temp_c'];
	this.enz_v = 0;
	this.buff_v = 0;
	this.dntp_v = 0;
	this.prim_t_v = 0;
	this.prim_b_v = 0;
	this.temp_v = 0;
	this.water_v = 0;
	this.total_v = 0;
}

MixFragment.prototype.go = function () {
    var m = this.mix.m();
    this.enz_v = m*this.mix.enz_d.value/this.mix.enz_s.value;
    var buff = m * this.mix.buff_d.value / this.mix.buff_s.value;
    var dntp = m * this.mix.dntp_d.value / this.mix.dntp_s.value;
    var primer_top = m * this.mix.prim_d.value / this.prim_t_s.value;
    var primer_bottom = m * this.mix.prim_d.value / this.prim_b_s.value;
    var template = m * this.mix.temp_d.value / this.temp_s.value;
    var total = (- template - this.enz_v) / (buff + dntp + primer_top + primer_bottom - 1);
    if (total < this.mix.vol_e.value) {
        this.total_v = this.mix.vol_e.value;
        this.water_v = this.total_v - this.total_v * (buff + dntp + primer_top + primer_bottom) - template - this.enz_v;
    } else {
        this.total_v = total;
        this.water_v = 0;
    }

	this.buff_v = this.total_v * buff;
	this.dntp_v = this.total_v * dntp;
	this.temp_v = template;
	this.prim_t_v = this.total_v * primer_top;
	this.prim_b_v = this.total_v * primer_bottom;
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
		if (this[x].constructor == Number){
			if (this[x] < 1 && this[x] > 0){
				//jQuery(this[x].parentNode).addClass("warningcell");
				this.mix.warnings += 1;
			} else {
				//jQuery(this[x].parentNode).removeClass("warningcell");
			}
		}
	}
}

jQuery('document').ready(function () {
	jQuery('#protocol-accordion').accordion({
		collapsible:true,
	});
	jQuery('#warning').hide();
    var mix = new Mix();
	jQuery('.mix').change(function(){
		validate(this);
		if (jQuery('.errorcell').size() > 0){
			jQuery('button#download_protocol').button('disable');
		} else {
            mix.go();
			jQuery('button#download_protocol').button('enable');
		}
	});
});
