# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'ConstructFragment.buffer_v'
        db.add_column(u'gibson_constructfragment', 'buffer_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.dntp_v'
        db.add_column(u'gibson_constructfragment', 'dntp_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.primer_top_v'
        db.add_column(u'gibson_constructfragment', 'primer_top_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.primer_bottom_v'
        db.add_column(u'gibson_constructfragment', 'primer_bottom_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.template_v'
        db.add_column(u'gibson_constructfragment', 'template_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.water_v'
        db.add_column(u'gibson_constructfragment', 'water_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)

        # Adding field 'ConstructFragment.total_v'
        db.add_column(u'gibson_constructfragment', 'total_v',
                      self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=3, decimal_places=1),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'ConstructFragment.buffer_v'
        db.delete_column(u'gibson_constructfragment', 'buffer_v')

        # Deleting field 'ConstructFragment.dntp_v'
        db.delete_column(u'gibson_constructfragment', 'dntp_v')

        # Deleting field 'ConstructFragment.primer_top_v'
        db.delete_column(u'gibson_constructfragment', 'primer_top_v')

        # Deleting field 'ConstructFragment.primer_bottom_v'
        db.delete_column(u'gibson_constructfragment', 'primer_bottom_v')

        # Deleting field 'ConstructFragment.template_v'
        db.delete_column(u'gibson_constructfragment', 'template_v')

        # Deleting field 'ConstructFragment.water_v'
        db.delete_column(u'gibson_constructfragment', 'water_v')

        # Deleting field 'ConstructFragment.total_v'
        db.delete_column(u'gibson_constructfragment', 'total_v')


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
            'buffer_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '100', 'max_digits': '4', 'decimal_places': '1'}),
            'construct': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': u"orm['gibson.Construct']"}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'dntp_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'end_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'end_feature'", 'null': 'True', 'to': u"orm['fragment.Feature']"}),
            'end_offset': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': u"orm['fragment.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'order': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'primer_bottom_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'primer_top_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'start_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'start_feature'", 'null': 'True', 'to': u"orm['fragment.Feature']"}),
            'start_offset': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'template_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'total_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'}),
            'water_v': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '3', 'decimal_places': '1'})
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

    complete_apps = ['gibson']