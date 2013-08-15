# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Cds'
        db.create_table(u'databaseInput_cds', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('origin', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['databaseInput.Origin'])),
            ('geneName', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('dnaSequence', self.gf('django.db.models.fields.TextField')()),
            ('description', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'databaseInput', ['Cds'])

        # Adding model 'Origin'
        db.create_table(u'databaseInput_origin', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sourceType', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('source', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('species', self.gf('django.db.models.fields.CharField')(max_length=100, null=True, blank=True)),
            ('product', self.gf('django.db.models.fields.CharField')(max_length=20, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'child', null=True, to=orm['databaseInput.Origin'])),
        ))
        db.send_create_signal(u'databaseInput', ['Origin'])

        # Adding model 'Domain'
        db.create_table(u'databaseInput_domain', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('module', self.gf('django.db.models.fields.IntegerField')()),
            ('cds', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['databaseInput.Cds'])),
            ('domainType', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['databaseInput.Type'])),
            ('chirality', self.gf('django.db.models.fields.CharField')(default=u'None', max_length=1)),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('pfamLinkerStart', self.gf('django.db.models.fields.IntegerField')()),
            ('pfamLinkerStop', self.gf('django.db.models.fields.IntegerField')()),
            ('definedLinkerStart', self.gf('django.db.models.fields.IntegerField')()),
            ('definedLinkerStop', self.gf('django.db.models.fields.IntegerField')()),
            ('pfamStart', self.gf('django.db.models.fields.IntegerField')()),
            ('pfamStop', self.gf('django.db.models.fields.IntegerField')()),
            ('definedStart', self.gf('django.db.models.fields.IntegerField')()),
            ('definedStop', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'databaseInput', ['Domain'])

        # Adding M2M table for field substrateSpecificity on 'Domain'
        m2m_table_name = db.shorten_name(u'databaseInput_domain_substrateSpecificity')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('domain', models.ForeignKey(orm[u'databaseInput.domain'], null=False)),
            ('substrate', models.ForeignKey(orm[u'databaseInput.substrate'], null=False))
        ))
        db.create_unique(m2m_table_name, ['domain_id', 'substrate_id'])

        # Adding model 'Substrate'
        db.create_table(u'databaseInput_substrate', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('chirality', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('structure', self.gf('django.db.models.fields.TextField')()),
            ('enantiomer', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['databaseInput.Substrate'], null=True, blank=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'child', null=True, to=orm['databaseInput.Substrate'])),
        ))
        db.send_create_signal(u'databaseInput', ['Substrate'])

        # Adding M2M table for field modification on 'Substrate'
        m2m_table_name = db.shorten_name(u'databaseInput_substrate_modification')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('substrate', models.ForeignKey(orm[u'databaseInput.substrate'], null=False)),
            ('modification', models.ForeignKey(orm[u'databaseInput.modification'], null=False))
        ))
        db.create_unique(m2m_table_name, ['substrate_id', 'modification_id'])

        # Adding model 'Modification'
        db.create_table(u'databaseInput_modification', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('domainType', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'modificationAdded', null=True, to=orm['databaseInput.Type'])),
        ))
        db.send_create_signal(u'databaseInput', ['Modification'])

        # Adding model 'Type'
        db.create_table(u'databaseInput_type', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('isModification', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('pfamName', self.gf('django.db.models.fields.CharField')(max_length=100, null=True, blank=True)),
            ('pfamId', self.gf('django.db.models.fields.CharField')(max_length=20, null=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'databaseInput', ['Type'])

        # Adding model 'Linkout'
        db.create_table(u'databaseInput_linkout', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('linkoutType', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['databaseInput.LinkoutType'])),
            ('identifier', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('content_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['contenttypes.ContentType'])),
            ('object_id', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'databaseInput', ['Linkout'])

        # Adding model 'LinkoutType'
        db.create_table(u'databaseInput_linkouttype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('shortcut', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('url', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('description', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'databaseInput', ['LinkoutType'])


    def backwards(self, orm):
        # Deleting model 'Cds'
        db.delete_table(u'databaseInput_cds')

        # Deleting model 'Origin'
        db.delete_table(u'databaseInput_origin')

        # Deleting model 'Domain'
        db.delete_table(u'databaseInput_domain')

        # Removing M2M table for field substrateSpecificity on 'Domain'
        db.delete_table(db.shorten_name(u'databaseInput_domain_substrateSpecificity'))

        # Deleting model 'Substrate'
        db.delete_table(u'databaseInput_substrate')

        # Removing M2M table for field modification on 'Substrate'
        db.delete_table(db.shorten_name(u'databaseInput_substrate_modification'))

        # Deleting model 'Modification'
        db.delete_table(u'databaseInput_modification')

        # Deleting model 'Type'
        db.delete_table(u'databaseInput_type')

        # Deleting model 'Linkout'
        db.delete_table(u'databaseInput_linkout')

        # Deleting model 'LinkoutType'
        db.delete_table(u'databaseInput_linkouttype')


    models = {
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
            'origin': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Origin']"})
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
            'substrateSpecificity': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['databaseInput.Substrate']", 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.linkout': {
            'Meta': {'object_name': 'Linkout'},
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'identifier': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'linkoutType': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.LinkoutType']"}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {})
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
            'product': ('django.db.models.fields.CharField', [], {'max_length': '20', 'null': 'True', 'blank': 'True'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'sourceType': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'species': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'})
        },
        u'databaseInput.substrate': {
            'Meta': {'object_name': 'Substrate'},
            'chirality': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'enantiomer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['databaseInput.Substrate']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modification': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['databaseInput.Modification']", 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'child'", 'null': 'True', 'to': u"orm['databaseInput.Substrate']"}),
            'structure': ('django.db.models.fields.TextField', [], {})
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