# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'Substrate.user'
        db.add_column(u'databaseInput_substrate', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)

        # Adding field 'Origin.user'
        db.add_column(u'databaseInput_origin', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)

        # Adding field 'Cds.user'
        db.add_column(u'databaseInput_cds', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)

        # Adding field 'Product.user'
        db.add_column(u'databaseInput_product', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)

        # Adding field 'Domain.user'
        db.add_column(u'databaseInput_domain', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)

        # Adding field 'Linkout.user'
        db.add_column(u'databaseInput_linkout', 'user',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'Substrate.user'
        db.delete_column(u'databaseInput_substrate', 'user_id')

        # Deleting field 'Origin.user'
        db.delete_column(u'databaseInput_origin', 'user_id')

        # Deleting field 'Cds.user'
        db.delete_column(u'databaseInput_cds', 'user_id')

        # Deleting field 'Product.user'
        db.delete_column(u'databaseInput_product', 'user_id')

        # Deleting field 'Domain.user'
        db.delete_column(u'databaseInput_domain', 'user_id')

        # Deleting field 'Linkout.user'
        db.delete_column(u'databaseInput_linkout', 'user_id')


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
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
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
            'product': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Product']", 'null': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.domain': {
            'Meta': {'object_name': 'Domain'},
            'cds': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Cds']"}),
            'chirality': ('django.db.models.fields.CharField', [], {'default': "u'None'", 'max_length': '1'}),
            'definedLinkerStart': ('django.db.models.fields.IntegerField', [], {}),
            'definedLinkerStop': ('django.db.models.fields.IntegerField', [], {}),
            'definedStart': ('django.db.models.fields.IntegerField', [], {}),
            'definedStop': ('django.db.models.fields.IntegerField', [], {}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'domainType': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Type']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'module': ('django.db.models.fields.IntegerField', [], {}),
            'pfamLinkerStart': ('django.db.models.fields.IntegerField', [], {}),
            'pfamLinkerStop': ('django.db.models.fields.IntegerField', [], {}),
            'pfamStart': ('django.db.models.fields.IntegerField', [], {}),
            'pfamStop': ('django.db.models.fields.IntegerField', [], {}),
            'substrateSpecificity': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['databaseInput.Substrate']", 'null': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.linkout': {
            'Meta': {'object_name': 'Linkout'},
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'identifier': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'linkoutType': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.LinkoutType']"}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.linkouttype': {
            'Meta': {'object_name': 'LinkoutType'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'shortcut': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'url': ('django.db.models.fields.CharField', [], {'max_length': '100'})
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
            'structure': ('django.db.models.fields.TextField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.type': {
            'Meta': {'object_name': 'Type'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isModification': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'pfamId': ('django.db.models.fields.CharField', [], {'max_length': '20', 'null': 'True', 'blank': 'True'}),
            'pfamName': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['databaseInput']