from django.contrib import admin
from django.contrib.contenttypes import generic

from databaseInput.models import (Domain, Substrate, 
	Origin, Type, Cds,	 Linkout, LinkoutType, Modification, Product,
	DomainTuple, Experiment)
# Register your models here.

class DomainInLine(generic.GenericTabularInline):
	model = Linkout
	extra = 1

class DomainAdmin(admin.ModelAdmin):
	list_display = ('domainType' , 'cds', 'module')
	inlines = [DomainInLine]


class CdsInLine(generic.GenericTabularInline):
	model = Linkout
	extra = 1

class ProductInLine(generic.GenericTabularInline):
    model = Linkout
    extra = 1

class ProductAdmin(admin.ModelAdmin):
    inlines = [ProductInLine]

class CdsAdmin(admin.ModelAdmin):
	inlines = [CdsInLine]

class OriginInLine(generic.GenericTabularInline):
	model = Linkout
	extra = 1

class OriginAdmin(admin.ModelAdmin):
	inlines = [OriginInLine]

class SubstrateInLine(generic.GenericTabularInline):
	model = Linkout
	extra = 1

class SubstrateAdmin(admin.ModelAdmin):
	inlines = [SubstrateInLine]

class ExperimentInLine(generic.GenericTabularInline):
	model = Linkout
	extra = 1

class ExperimentAdmin(admin.ModelAdmin):
	inlines = [SubstrateInLine]

admin.site.register(Domain, DomainAdmin)
admin.site.register(Substrate, SubstrateAdmin)
admin.site.register(Origin, OriginAdmin)
admin.site.register(Cds, CdsAdmin)
admin.site.register(Type)
admin.site.register(Product, ProductAdmin)
admin.site.register(LinkoutType)
admin.site.register(Linkout)
admin.site.register(Modification)
admin.site.register(Experiment,ExperimentAdmin)
admin.site.register(DomainTuple)

