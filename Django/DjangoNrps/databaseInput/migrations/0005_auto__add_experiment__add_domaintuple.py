# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Experiment'
        db.create_table(u'databaseInput_experiment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True)),
        ))
        db.send_create_signal(u'databaseInput', ['Experiment'])

        # Adding model 'DomainTuple'
        db.create_table(u'databaseInput_domaintuple', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('domain_left', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'next_domains', to=orm['databaseInput.Domain'])),
            ('domain_right', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'previous_domains', to=orm['databaseInput.Domain'])),
            ('position_left', self.gf('django.db.models.fields.IntegerField')()),
            ('position_right', self.gf('django.db.models.fields.IntegerField')()),
            ('experiment', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'domain_tuples', to=orm['databaseInput.Experiment'])),
        ))
        db.send_create_signal(u'databaseInput', ['DomainTuple'])


    def backwards(self, orm):
        # Deleting model 'Experiment'
        db.delete_table(u'databaseInput_experiment')

        # Deleting model 'DomainTuple'
        db.delete_table(u'databaseInput_domaintuple')


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
            'product': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Product']", 'null': 'True', 'blank': 'True'}),
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
            'pfamLinkerStart': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pfamLinkerStop': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pfamStart': ('django.db.models.fields.IntegerField', [], {}),
            'pfamStop': ('django.db.models.fields.IntegerField', [], {}),
            'substrateSpecificity': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "u'adenylationDomain'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['databaseInput.Substrate']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.domaintuple': {
            'Meta': {'object_name': 'DomainTuple'},
            'domain_left': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'next_domains'", 'to': u"orm['databaseInput.Domain']"}),
            'domain_right': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'previous_domains'", 'to': u"orm['databaseInput.Domain']"}),
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'domain_tuples'", 'to': u"orm['databaseInput.Experiment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'position_left': ('django.db.models.fields.IntegerField', [], {}),
            'position_right': ('django.db.models.fields.IntegerField', [], {})
        },
        u'databaseInput.experiment': {
            'Meta': {'object_name': 'Experiment'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
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
        }
    }

    complete_apps = ['databaseInput']