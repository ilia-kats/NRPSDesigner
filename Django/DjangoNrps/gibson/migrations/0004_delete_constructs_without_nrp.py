# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import DataMigration
from django.db import models
from django.core.exceptions import ObjectDoesNotExist

class Migration(DataMigration):

    def forwards(self, orm):
        "Write your forwards methods here."
        # Note: Remember to use orm['appname.ModelName'] rather than "from appname.models..."
        for c in orm.Construct.objects.filter(nrp__isnull=True):
            for x in c.cf.all():
                if x.fragment is not None and x.fragment.origin == 'ND':
                    try:
                        g = orm['fragment.Gene'].objects.get(owner=x.fragment.owner, pk=x.fragment.pk)
                        try:
                            orm['fragment.Reference'].objects.filter(gene=g).delete()
                        except ObjectDoesNotExist:
                            pass
                        try:
                            orm['fragment.Annotation'].objects.filter(gene=g).delete()
                        except ObjectDoesNotExist:
                            pass
                        try:
                            feats = orm['fragment.Feature'].objects.filter(gene=g)
                            for f in feats:
                                try:
                                    orm['fragment.Qualifier'].objects.filter(feature=f).delete()
                                except ObjectDoesNotExist:
                                    pass
                            feats.delete()
                        except ObjectDoesNotExist:
                            pass
                        g.delete()
                    except ObjectDoesNotExist:
                        pass
                x.delete()
            if c.primer is not None:
                [x.del_all() for x in c.primer.all() if x is not None]
            c.delete()

    def backwards(self, orm):
        "Write your backwards methods here."

    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'databaseInput.cds': {
            'Meta': {'object_name': 'Cds'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'dnaSequence': ('django.db.models.fields.TextField', [], {}),
            'geneName': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'origin': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Origin']"}),
            'product': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'cdsSequence'", 'null': 'True', 'to': u"orm['databaseInput.Product']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.domain': {
            'Meta': {'object_name': 'Domain'},
            'cds': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'domains'", 'to': u"orm['databaseInput.Cds']"}),
            'chirality': ('django.db.models.fields.CharField', [], {'default': "u'N'", 'max_length': '1'}),
            'definedLinkerStart': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'definedLinkerStop': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'definedStart': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'definedStop': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'domainType': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Type']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'module': ('django.db.models.fields.IntegerField', [], {}),
            'next_domain': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'prev_domain'", 'blank': 'True', 'through': u"orm['databaseInput.DomainTuple']", 'to': u"orm['databaseInput.Domain']"}),
            'pfamLinkerStart': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pfamLinkerStop': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pfamStart': ('django.db.models.fields.IntegerField', [], {}),
            'pfamStop': ('django.db.models.fields.IntegerField', [], {}),
            'substrateSpecificity': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "u'adenylationDomain'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['databaseInput.Substrate']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.domaintuple': {
            'Meta': {'object_name': 'DomainTuple'},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'domain_tuples'", 'to': u"orm['databaseInput.Experiment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'next_domain': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'prev_tuple'", 'to': u"orm['databaseInput.Domain']"}),
            'next_position': ('django.db.models.fields.IntegerField', [], {}),
            'prev_domain': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'next_tuple'", 'to': u"orm['databaseInput.Domain']"}),
            'prev_position': ('django.db.models.fields.IntegerField', [], {})
        },
        u'databaseInput.experiment': {
            'Meta': {'object_name': 'Experiment'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.modification': {
            'Meta': {'object_name': 'Modification'},
            'domainType': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'modificationAdded'", 'null': 'True', 'to': u"orm['databaseInput.Type']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'databaseInput.origin': {
            'Meta': {'object_name': 'Origin'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'child'", 'null': 'True', 'to': u"orm['databaseInput.Origin']"}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'sourceType': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.product': {
            'Meta': {'object_name': 'Product'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.substrate': {
            'Meta': {'object_name': 'Substrate'},
            'chirality': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'enantiomer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Substrate']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['databaseInput.Modification']", 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'child'", 'null': 'True', 'to': u"orm['databaseInput.Substrate']"}),
            'smashName': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'structure': ('django.db.models.fields.TextField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.type': {
            'Meta': {'object_name': 'Type'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isModification': ('django.db.models.fields.BooleanField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'pfamGraphic': ('django.db.models.fields.TextField', [], {}),
            'pfamId': ('django.db.models.fields.CharField', [], {'max_length': '20', 'null': 'True', 'blank': 'True'}),
            'pfamName': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'smashName': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'})
        },
        u'designerGui.domainorder': {
            'Meta': {'object_name': 'DomainOrder'},
            'domain': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'domainOrder'", 'to': u"orm['databaseInput.Domain']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'nrp': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'domainOrder'", 'to': u"orm['designerGui.NRP']"}),
            'order': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'designerGui.nrp': {
            'Meta': {'object_name': 'NRP'},
            'construct': ('django.db.models.fields.related.OneToOneField', [], {'blank': 'True', 'related_name': "'nrp'", 'unique': 'True', 'null': 'True', 'to': u"orm['gibson.Construct']"}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '2000', 'null': 'True', 'blank': 'True'}),
            'designed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'designerDomains': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "'includedIn'", 'blank': 'True', 'through': u"orm['designerGui.DomainOrder']", 'to': u"orm['databaseInput.Domain']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'indigoidineTagged': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'monomers': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "'includedIn'", 'blank': 'True', 'through': u"orm['designerGui.SubstrateOrder']", 'to': u"orm['databaseInput.Substrate']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'child'", 'null': 'True', 'to': u"orm['designerGui.NRP']"}),
            'uuid': ('django.db.models.fields.CharField', [], {'default': "'6d1b367c-4ce3-4b5f-9386-0cc5f8b455d9'", 'max_length': '36', 'db_index': 'True'})
        },
        u'designerGui.species': {
            'Meta': {'object_name': 'Species'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'taxon_id': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        u'designerGui.substrateorder': {
            'Meta': {'object_name': 'SubstrateOrder'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'nrp': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'substrateOrder'", 'to': u"orm['designerGui.NRP']"}),
            'order': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'substrate': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'substrateOrder'", 'to': u"orm['databaseInput.Substrate']"})
        },
        u'fragment.annotation': {
            'Meta': {'object_name': 'Annotation'},
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'annotations'", 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '5120', 'blank': 'True'})
        },
        u'fragment.domaingene': {
            'Meta': {'object_name': 'DomainGene', '_ormbases': [u'fragment.Gene']},
            'domains': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'+'", 'symmetrical': 'False', 'to': u"orm['databaseInput.Domain']"}),
            u'gene_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['fragment.Gene']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'fragment.feature': {
            'Meta': {'ordering': "['start']", 'object_name': 'Feature'},
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'end': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'features'", 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'start': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        u'fragment.gene': {
            'Meta': {'object_name': 'Gene'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'origin': ('django.db.models.fields.CharField', [], {'max_length': '2'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'max_length': '500000'}),
            'viewable': ('django.db.models.fields.CharField', [], {'default': "'L'", 'max_length': '1'})
        },
        u'fragment.qualifier': {
            'Meta': {'object_name': 'Qualifier'},
            'data': ('django.db.models.fields.CharField', [], {'max_length': '512'}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'qualifiers'", 'to': u"orm['fragment.Feature']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        u'fragment.reference': {
            'Meta': {'object_name': 'Reference'},
            'authors': ('django.db.models.fields.CharField', [], {'max_length': '1024'}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'references'", 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'journal': ('django.db.models.fields.CharField', [], {'max_length': '512'}),
            'medline_id': ('django.db.models.fields.CharField', [], {'max_length': '24', 'blank': 'True'}),
            'pubmed_id': ('django.db.models.fields.CharField', [], {'max_length': '24', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '1024'})
        },
        u'gibson.construct': {
            'Meta': {'object_name': 'Construct'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '2000'}),
            'fragments': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "'construct_slave'", 'blank': 'True', 'through': u"orm['gibson.ConstructFragment']", 'to': u"orm['fragment.Gene']"}),
            'genbank': ('django.db.models.fields.related.OneToOneField', [], {'blank': 'True', 'related_name': "'construct_master'", 'unique': 'True', 'null': 'True', 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True'}),
            'processed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'shape': ('django.db.models.fields.CharField', [], {'max_length': '1'})
        },
        u'gibson.constructfragment': {
            'Meta': {'ordering': "['order']", 'object_name': 'ConstructFragment'},
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '100', 'max_digits': '4', 'decimal_places': '1'}),
            'construct': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': u"orm['gibson.Construct']"}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'end_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'end_feature'", 'null': 'True', 'to': u"orm['fragment.Feature']"}),
            'end_offset': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'order': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'start_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'start_feature'", 'null': 'True', 'to': u"orm['fragment.Feature']"}),
            'start_offset': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        },
        u'gibson.pcrsettings': {
            'Meta': {'object_name': 'PCRSettings'},
            'buffer_d': ('django.db.models.fields.DecimalField', [], {'default': '1', 'max_digits': '3', 'decimal_places': '1'}),
            'buffer_s': ('django.db.models.fields.DecimalField', [], {'default': '10', 'max_digits': '3', 'decimal_places': '1'}),
            'construct': ('annoying.fields.AutoOneToOneField', [], {'related_name': "'pcrsettings'", 'unique': 'True', 'to': u"orm['gibson.Construct']"}),
            'dntp_d': ('django.db.models.fields.DecimalField', [], {'default': '0.8', 'max_digits': '3', 'decimal_places': '1'}),
            'dntp_s': ('django.db.models.fields.DecimalField', [], {'default': '10', 'max_digits': '3', 'decimal_places': '1'}),
            'enzyme_d': ('django.db.models.fields.DecimalField', [], {'default': '2.5', 'max_digits': '3', 'decimal_places': '1'}),
            'enzyme_s': ('django.db.models.fields.DecimalField', [], {'default': '2.5', 'max_digits': '3', 'decimal_places': '1'}),
            'error_margin': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '10'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'primer_d': ('django.db.models.fields.DecimalField', [], {'default': '0.4', 'max_digits': '3', 'decimal_places': '1'}),
            'repeats': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '1'}),
            'template_d': ('django.db.models.fields.DecimalField', [], {'default': '100', 'max_digits': '4', 'decimal_places': '1'}),
            'volume_each': ('django.db.models.fields.DecimalField', [], {'default': '12.5', 'max_digits': '3', 'decimal_places': '1'})
        },
        u'gibson.primer': {
            'Meta': {'ordering': "['stick']", 'object_name': 'Primer'},
            'boxplot': ('django.db.models.fields.files.ImageField', [], {'max_length': '100'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '5', 'max_digits': '4', 'decimal_places': '1'}),
            'construct': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'primer'", 'to': u"orm['gibson.Construct']"}),
            'flap': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'flap'", 'unique': 'True', 'to': u"orm['gibson.PrimerHalf']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'stick': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'stick'", 'unique': 'True', 'to': u"orm['gibson.PrimerHalf']"})
        },
        u'gibson.primerhalf': {
            'Meta': {'object_name': 'PrimerHalf'},
            'cfragment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'ph'", 'to': u"orm['gibson.ConstructFragment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'top': ('django.db.models.fields.BooleanField', [], {})
        },
        u'gibson.settings': {
            'Meta': {'object_name': 'Settings'},
            'construct': ('annoying.fields.AutoOneToOneField', [], {'related_name': "'settings'", 'unique': 'True', 'to': u"orm['gibson.Construct']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mg_salt': ('django.db.models.fields.DecimalField', [], {'default': '0.0', 'max_digits': '3', 'decimal_places': '2'}),
            'min_anneal_tm': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '50'}),
            'min_overlap': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '20'}),
            'min_primer_tm': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '60'}),
            'na_salt': ('django.db.models.fields.DecimalField', [], {'default': '0.05', 'max_digits': '4', 'decimal_places': '3'}),
            'ss_safety': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '3'})
        },
        u'gibson.warning': {
            'Meta': {'object_name': 'Warning'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'primer': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'warning'", 'to': u"orm['gibson.Primer']"}),
            'text': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '2'})
        }
    }

    complete_apps = ['designerGui', 'fragment', 'gibson']
    symmetrical = True
