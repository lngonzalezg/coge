define([ 'dojo/_base/declare',
         'dojo/_base/array',
         'JBrowse/View/Export'
       ],
       function( declare, array, ExportBase ) {

return declare( ExportBase,
 /**
  * @lends JBrowse.View.Export.GFF3.prototype
  */
{

    /**
     * Data export driver for GFF3 format.
     * @constructs
     */
    constructor: function( args ) {
        this._idCounter = 0;

        this.print( "##gff-version 3\n");
        if( this.refSeq )
            this.print( "##sequence-region "+this.refSeq.name+" "+(this.refSeq.start+1)+" "+this.refSeq.end+"\n" );

        this.lastSync = 0;
    },

    gff3_field_names: [
        'seq_id',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'attributes'
    ],

    gff3_reserved_attributes: [
        'ID',
        'Name',
        'Alias',
        'Parent',
        'Target',
        'Gap',
        'Derives_from',
        'Note',
        'Dbxref',
        'Ontology_term',
        'Is_circular'
    ],

    /**
     * @returns false if the field goes in tabular portion of gff3, true otherwise
     * @private
     */
    _is_not_gff3_tab_field: function( fieldname ) {
        if( ! this._gff3_fields_by_name ) {
            var fields = {};
            dojo.forEach( this.gff3_field_names, function(f) {
                              fields[f] = true;
                          });
            this._gff3_fields_by_name = fields;
        }

        return ! this._gff3_fields_by_name[ fieldname.toLowerCase() ];
    },

    /**
     * @returns the capitalized attribute name if the given field name
     * corresponds to a GFF3 reserved attribute
     * @private
     */
    _gff3_reserved_attribute: function( fieldname ) {
        if( ! this._gff3_reserved_attributes_by_lcname ) {
            var fields = {};
            dojo.forEach( this.gff3_reserved_attributes, function(f) {
                              fields[f.toLowerCase()] = f;
                          });
            this._gff3_reserved_attributes_by_lcname = fields;
        }

        return this._gff3_reserved_attributes_by_lcname[ fieldname.toLowerCase() ];
    },


    /**
     * Format a feature into a string.
     * @param {Object} feature feature object (like those returned from JBrowse/Store/SeqFeature/*)
     * @returns {String} GFF3 string representation of the feature
     */
    formatFeature: function( feature, parentID ) {
        var fields = dojo.map(
                [ feature.get('seq_id') || this.refSeq.name ]
                .concat( dojo.map( this.gff3_field_names.slice(1,7), function(field) {
                                       return feature.get( field );
                                   },this)
                       ),
            function( data ) {
                var dt = typeof data;
                return this._gff3_escape( dt == 'string' || dt == 'number' ? data : '.' );
            },
            this
        );

        // convert back from interbase
        if( typeof parseInt(fields[3]) == 'number' )
            fields[3]++;
        // normalize the strand field
        fields[6] = { '1': '+', '-1': '-', '0': '.' }[ fields[6] ] || fields[6];

        // format the attributes
        var attr = this._gff3_attributes( feature );
        if( parentID )
            attr.Parent = parentID;
        else
            delete attr.Parent;

        var subfeatures = array.map(
            feature.get('subfeatures') || [],
            function(feat) {
                if( ! attr.ID ) {
                    attr.ID = ++this._idCounter;
                }
                return this.formatFeature( feat, attr.ID );
            }, this);

        // need to format the attrs after doing the subfeatures,
        // because the subfeature formatting might have autocreated an
        // ID for the parent
        fields[8] = this._gff3_format_attributes( attr );

        return [ fields.join("\t")+"\n" ].concat( subfeatures );
    },

    /**
     * Write the feature to the GFF3 under construction.
     * @returns nothing
     */
    writeFeature: function(feature) {
        var fmt = this.formatFeature(feature);
        this.print( fmt );

        // avoid printing sync marks more than every 10 lines
        if( this.lastSync >= 9 ) {
            this.lastSync = 0;
            this.print( "###\n" );
        } else {
            this.lastSync += fmt.length || 1;
        }
    },

    /**
     * Extract a key-value object of gff3 attributes from the given
     * feature.  Attribute names will have proper capitalization.
     * @private
     */
    _gff3_attributes: function(feature) {
        var tags = array.filter( feature.tags(), dojo.hitch(this, function(f) {
            f = f.toLowerCase();
            return this._is_not_gff3_tab_field(f) && f != 'subfeatures';
        }));
        var attrs = {};
        array.forEach( tags, function(tag) {
            var val = feature.get(tag);
            var valtype = typeof val;
            if( valtype == 'boolean' )
                val = val ? 1 : 0;
            else if( valtype == 'undefined' )
                return;
            tag = this._gff3_reserved_attribute(tag) || this._ensure_non_reserved( tag );
            attrs[tag] = val;
        },this);
        return attrs;
    },

    // ensure that an attribute name is not reserved.  currently does
    // this by adding a leading underscore to attribute names that
    // have initial capital letters.
    _ensure_non_reserved: function( str ) {
        return str.replace(/^[A-Z]/,function() { return '_'+str[0]; });
    },

    /**
     * @private
     * @returns {String} formatted attribute string
     */
    _gff3_format_attributes: function( attrs ) {
        var attrOrder = [];
        for( var tag in attrs ) {
            attrOrder.push( this._gff3_escape( tag )+'='+this._gff3_escape( attrs[tag] ) );
        }
        return attrOrder.join(';');
    },

    /**
     * @returns always an escaped string representation of the passed value
     * @private
     */
    _gff3_escape: function( val ) {
        return (''+val).replace(/[\n\r\t\;\=%&,\x00-\x1f\x7f-\xff]+/g, escape );
    }
});

});
