/*
 * Fragment API
 *  --fetch fragment from the server
 *
 */

var libFrag = new function()
{
	/*
	 * Public functions
	 *
	 */

	//fetch a fragment by its ID
	this.getByID = function(id, _suc)
	{
		AJAX.post({
			url: fragment_api_url + id+'/',
			success: function(data)
			{
				_suc(new Fragment(data));
			},
		});
	};

    //Fetch all available fragments, returned in a list
	this.getAll = function()
	{
        var _suc = Array.prototype.slice.call(arguments,-1)[0];
        postdata = {
            url: fragment_api_url + 'listAll/',
            success: function(data)
            {
                var frags = new Array();
                for(f in data)
                {
                    frags.push(new Fragment(data[f]));
                }
                _suc(frags);
            },
        };
        if (arguments.length > 1)
            postdata['data'] = {'type': Array.prototype.slice.call(arguments, 0, -1)};
        AJAX.post(postdata);
	}

    /*
     * Return fragment colors
     */
    var _h = Math.random() * 360;
	var _grc = 0.618033988749895 * 360;
	
	this.getNextColor = function()
	{
		_h = (_h + _grc) % 360;
		return Graphics.getHSL(_h,40,50);
	}


	//IUAPC Unambiguous Alphabet / Complement table
	this.alphabet = {
		'A': 'T',
		'B': 'V',
		'C': 'G',
		'D': 'H',
		'G': 'C',
		'H': 'D',
		'K': 'M',
		'M': 'K',
		'N': 'N',
		'R': 'Y',
		'S': 'S',
		'T': 'A',
		'V': 'B',
		'W': 'W',
		'Y': 'R',
		'a': 't',
		'b': 'v',
		'c': 'g',
		'd': 'h',
		'g': 'c',
		'h': 'd',
		'k': 'm',
		'm': 'k',
		'n': 'n',
		'r': 'y',
		's': 's',
		't': 'a',
		'v': 'b',
		'w': 'w',
		'y': 'r',
		' ': ' ',
	}
}

/*
 *
 * An actual fragment object
 *
 */
function Fragment(data)
{
	/*
	 * Public Accessors
	 *
	 */
	this.getID = function() {return id;};
	this.getName = function() {return name;};
    this.getDesc = function() {return desc;};
	this.getLength = function() {return length;};

	//if the sequence has not already been fetched, it will be streamed, otherwise
	//complete_fn will be called immediately
	this.getSequence = function( update_fn, complete_fn)
	{
		if(sequence!=null)
		{
			complete_fn(sequence);
			return;
		}

		//stream the sequence from the server
		AJAX.stream({
			url: fragment_api_url + id + '/getSeq/',
			success: function(data)
			{
				sequence = new String();
				for(i in data)
				{
					if(data[i] != ' ')
						sequence = sequence + data[i];
				}
				complete_fn(sequence);
			},
			type: 'GET',
			error: function(error)
			{
				console.error('Error fetching sequence: ' + error);
			}
		},
		function(data)
		{
			sequence = new String();
			for(i in data)
				{
					if(data[i] != ' ')
						sequence = sequence + data[i];
				}
			update_fn(sequence);
		} );
	};

	//if the metadata has not already been fetched, it will be fetched otherwise
	//complete_fn will be called immediately
	this.getMeta = function(complete_fn)
	{
		if(metadata!=null)
		{
			complete_fn(metadata);
			return;
		}

		//fetch the metadata from the server
		AJAX.post({
			url: fragment_api_url + id + '/getMeta/',
			success: function(ret){
				metadata = ret;
				complete_fn(ret);
			},
			error: function() {},
		});
	};

	this.setMeta = function(new_meta, success_cb, fail_cb)
	{
		metadata = new_meta;
		AJAX.post({
			url: fragment_api_url + id + '/setMeta/',
			success: success_cb,
			error: fail_cb,
			data: metadata,
		});
	};	

	this.getFeats = function(success_fn)
	{
		AJAX.post({
			url: fragment_api_url + id + '/getFeats/',
			success: success_fn,
			error: function(err)
			{
				console.error('Error getting features for ' + name);
			},
		});
	};

    /*
     * cropsettings = {
     *  start: start position,
     *  end: end position,
     *  f_internal: keep internal features?,
     *  f_all: keep all features
     *  result: in ['new', 'overwrite']
     *  new_name: defined if result == new
     *  new_desc: defined if result == new
     *  };
     *
     *  */
    this.crop = function(cropsettings, cb)
    {
        AJAX.post({
            url: fragment_api_url + id + '/crop/',
            data: cropsettings,
            success: function(data){
                if(jQuery.isFunction(cb)) cb(data);
            },
            error: function() {
                console.error('Error cropping fragment');
            },
        });
    };

    this.toString = function()
    {
        return '[Fragment name="' + name + '"]';
    }

	/*
	 *
	 * Private data
	 *
	 */

	var id = data.id;
	var name = data.name;
	var desc = data.desc;
	var length = data.length;

	//Data which might get filled in by subsequent AJAX calls
	var sequence = null;
	var metadata = ('annots' in data && 'refs' in data) ? data : null;
}

// ============================================ Widgets

/*
 *  jFragment: draggable fragment
 *
 * */
jQuery.widget("ui.jFragment", jQuery.ui.draggable, {
    options: {
        fragment: null,
        color: 'red',
        draggable: true,
    },
    _create: function() {
        this.f = this.options.fragment;
        this.jQueryel = jQuery(this.element[0])
        if(this.options.draggable)
            this.jQueryel.draggable(this.options);
        this.jQueryel.html("<p class='jf-name'>"+this.f.getName()+"</p>");
        this.jQueryel.addClass('jFragment');
        this.jQueryel.css({'background-color':this.options.color,
                     'border-color':this.options.color,});
    },
    getFragment: function()
    {
        return this.f;
    },
    getName: function()
    {
        return this.f.getName();
    },
    getColor: function()
    {
        return this.jQueryel.css('background-color');
    },
    option: function(name, value)
    {
        if(value == undefined)
            return this.options[name];
        console.log('setting '+name);
        this.options[name] = value;
    },
});  

/*
 *
 * jFragmentSelector: Easily select fragment from the available ones 
 * using drag and drop
 *
 * */
jQuery.widget("ui.jFragmentSelector", {
    options: {
        droptarget: null,
        containment: 'parent',
        dragEnabled: true,
    },
    _create: function() {
        var self = this;
        var jQuerye = this.jQueryel = jQuery(this.element[0]);
        this.jQueryfragView = jQuerye.find('.JFS_fragView');
        this.jQueryfilterHolder = jQuerye.find('.JFS_filterHolder');
        this.jQueryfilterHint = this.jQueryfilterHolder.find('p');
        this.jQueryfilterInput = this.jQueryfilterHolder.find('input');
        this.jQuerynumItems = jQuerye.find('#JFS_num_items');

        this.filter_timeout = null;
        //if we click the hint, we select the input
        this.jQueryfilterHint.on('click', function() {
            self.jQueryfilterInput.focus();
        });
        this.jQueryfilterInput.on('focus', function() {
            self._onInputFocus();
        });
        this.jQueryfilterInput.on('blur', function() {
            self._onInputBlur();
        });
        this.jQueryfilterInput.keypress( function(){
            if(self.filter_timeout != null){
                clearTimeout(self.filter_timeout);
                self.filter_timeout = null;
            }
            self.filter_timeout = setTimeout(function(){
                self._filter();
            }, 350);
        });
        this.jQueryfilterInput.hover(function() {
            self.jQueryfilterHolder.addClass('ui-state-hover');   
        }, function() {
            self.jQueryfilterHolder.removeClass('ui-state-hover');   
        });

        this.frags = new Array();
        //fetch the fragments
        libFrag.getAll("Expression_Plasmid", "Plasmid_Backbone", "Promoter", "RBS", "Terminator", function(frags){
            var categories = new Array();
            for (f in frags)
            {
                frags[f].getMeta(function(meta){
                    if ('annots' in meta && 'part_type' in meta['annots'])
                    {
                        if (!(meta['annots']['part_type'] in categories))
                            categories[meta['annots']['part_type']] = new Array();
                        categories[meta['annots']['part_type']].push(frags[f]);
                    }
                })
            }
            //remove the loading screen
            self.jQueryfragView.empty();
            //add in the fragments one by one
            for (c in categories)
            {
                var header = jQuery('<h5>' + c + '</h5>');
                var div = jQuery('<div/>');
                for (f in categories[c])
                {
                    self.frags.push(categories[c][f]);
                    jQuery('<div/>').addClass('JFS_fragHolder')
                    .append( jQuery('<div/>').jFragment({
                        fragment: categories[c][f],
                        color: libFrag.getNextColor(),
                        helper: function(){
                            return jQuery('<div/>').jFragment({
                                draggable:false,
                                fragment: jQuery(this).jFragment('getFragment'),
                            }).appendTo(self.jQueryel);
                        },
                        containment: self.options.containment,
                        zIndex:200,
                    }))
                    .data('f', categories[c][f])
                    .appendTo(div);
                }
                header.appendTo(self.jQueryfragView);
                div.appendTo(self.jQueryfragView);
            }
            self.jQueryfragView.accordion({collapsible: true, heightStyle: "content"});
            self.jQuerynumItems.text(frags.length);
        });


    },
    //set height to fill the parent container
    height: function(h){
        this.jQueryfragView.height( h - 
                              (this.jQueryel.height() - this.jQueryfragView.height()));
    },
    _onInputFocus: function(){
        this.jQueryfilterHint.hide();
        this.jQueryfilterHolder.addClass('ui-state-focus');
    },
    _onInputBlur: function(){
        if(!this.jQueryfilterInput.val())
            this.jQueryfilterHint.show();
        this.jQueryfilterHolder.removeClass('ui-state-focus');
    },
    _filter: function(){
        var s = this.jQueryfilterInput.val();
        var self = this;
        var matches = 0;
        if(s){
            this.jQueryel.find('.JFS_fragHolder').each(function(index, el){
                //hide the element if it doesn't match s
                if( self.frags[index].getName().toLowerCase()
                   .indexOf(s.toLowerCase()) < 0)
                {
                    jQuery(this).hide();
                }
                else
                {
                    jQuery(this).show();
                    matches = matches + 1;
                }
            });
        }
        else{
            this.jQueryel.find('.JFS_fragHolder').each(function(index, el){
                //show all if s is empty
                jQuery(this).show();
                matches = matches + 1;
            });
        }
        this.jQuerynumItems.text(matches);
    },
    option: function(name, value)
    {
        if(value == undefined)
            return this.options[name];
        console.log('setting '+name);
        this.options[name] = value;
    },
});


