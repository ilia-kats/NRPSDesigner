from django.db import models

class Species(models.Model):
    species = models.CharField(max_length=100)
    taxon_id = models.CharField(max_length=20)
  
    def __unicode__(self):
        return str(species)

    class Meta:
        verbose_name = "Species"
        verbose_name_plural = "Species"
