# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'DomainGene'
        db.create_table(u'fragment_domaingene', (
            (u'gene_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['fragment.Gene'], unique=True, primary_key=True)),
        ))
        db.send_create_signal(u'fragment', ['DomainGene'])

        # Adding M2M table for field domains on 'DomainGene'
        db.create_table(u'fragment_domaingene_domains', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('domaingene', models.ForeignKey(orm[u'fragment.domaingene'], null=False)),
            ('domain', models.ForeignKey(orm[u'databaseInput.domain'], null=False))
        ))
        db.create_unique(u'fragment_domaingene_domains', ['domaingene_id', 'domain_id'])


    def backwards(self, orm):
        # Deleting model 'DomainGene'
        db.delete_table(u'fragment_domaingene')

        # Removing M2M table for field domains on 'DomainGene'
        db.delete_table('fragment_domaingene_domains')


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
        }
    }

    complete_apps = ['fragment']