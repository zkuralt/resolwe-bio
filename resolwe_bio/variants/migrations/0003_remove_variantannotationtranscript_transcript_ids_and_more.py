# Generated by Django 4.2.13 on 2024-09-02 13:13

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        (
            "resolwe_bio_variants",
            "0002_rename_unfiltered_allele_depth_variantcall_alternative_allele_depth",
        ),
    ]

    operations = [
        migrations.RemoveField(
            model_name="variantannotationtranscript",
            name="transcript_ids",
        ),
        migrations.AddField(
            model_name="variantannotationtranscript",
            name="transcript_id",
            field=models.CharField(default="", max_length=200),
            preserve_default=False,
        ),
    ]
