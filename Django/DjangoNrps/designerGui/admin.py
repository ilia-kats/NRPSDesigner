from django.contrib import admin
from designerGui.models import NRP, DomainOrder

class NRPInline(admin.TabularInline):
    model = NRP.monomers.through
    extra = 3

class NRPAdmin(admin.ModelAdmin):
	list_display = ('name' , 'owner', 'modified')
	inlines = [NRPInline]

admin.site.register(NRP, NRPAdmin)
admin.site.register(DomainOrder)
