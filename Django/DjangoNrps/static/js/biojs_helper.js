function BiojsSequenceHelper(options)
{
    this.opts = {
        columns: {size:100},
        format: "FASTA",
        formatSelectorVisible: false,
    };
    for (var opt in options) {
        this.opts[opt] = options[opt];
    }
    this.seq = new Biojs.Sequence(this.opts);

    if ('target' in this.opts) {
        var self = this;
        this.seq._container.ready(function(){
            var chr = jQuery("<span>&nbsp;</span>");
            var width = chr.appendTo(self.seq._contentDiv).outerWidth(true);
            self.opts.columns.size = Math.floor(self.seq._contentDiv.innerWidth() / width * 0.9);
            self.seq._contentDiv.hide();
            chr.remove();
        });
    }

    this.showBioJsSequence = function (data, scrollTo) {
        if (typeof scrollTo == 'undefined')
            scrollTo = this.scrollTo.start
        this.seq.clearSequence();
        var self = this;
        jQuery.get(this.opts.url,
            data,
            function(data) {
                if (data[0] !== -1) {
                    if (typeof(data[1].biojs) !== 'undefined') {
                        var options = self.opts;
                        for (var attr in data[1].biojs) {
                            options[attr] = data[1].biojs[attr];
                        }
                        self.seq = new Biojs.Sequence(options);
                        if (typeof(self.opts.onAnnotationClicked) !== 'undefined')
                            self.seq.onAnnotationClicked(self.opts.onAnnotationClicked);
                        if (typeof(self.opts.onSelectionChange) !== 'undefined')
                            self.seq.onSelectionChange(self.opts.onSelectionChange);
                        if (typeof(self.opts.onSelectionChanged) !== 'undefined')
                            self.seq.onSelectionChanged(self.opts.onSelectionChanged);
                    }
                    if (typeof(data[1].highlight) != 'undefined') {
                        self.seq.addHighlight(data[1].highlight);
                        var seqelems = self.seq._contentDiv.find('.sequence');
                        if (scrollTo == self.scrollTo.start)
                            self.seq._container.scrollTo(seqelems[data[1].highlight.start - 1]);
                        else if (scrollTo == self.scrollTo.end)
                            self.seq._container.scrollTo(seqelems[data[1].highlight.end]);
                    }
                }
        });
    };

    this.showBioJsDomain = function (domain, scrollTo) {
        if (typeof scrollTo == 'undefined')
            scrollTo = this.scrollTo.start
        this.seq.unHighlightAll();
        var self = this;
        jQuery.get(this.opts.url,
            {domainId: domain, domainOnly: true},
            function(data) {
                if (data[0] !== -1) {
                    if (typeof(data[1].highlight) != 'undefined') {
                        self.seq.addHighlight(data[1].highlight);
                        var seqelems = self.seq._contentDiv.find('.sequence');
                        if (scrollTo == self.scrollTo.start)
                            self.seq._container.scrollTo(seqelems[data[1].highlight.start - 1]);
                        else if (scrollTo == self.scrollTo.end)
                            self.seq._container.scrollTo(seqelems[data[1].highlight.end]);
                    }
                }
        });
    }

    this.getBiojsSequence = function() {
        return this.seq;
    }
}

BiojsSequenceHelper.prototype.scrollTo = {start: 0, end: 1, none: 2};
