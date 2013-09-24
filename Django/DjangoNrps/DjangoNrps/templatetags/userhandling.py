from django import template

register = template.Library()

@register.filter
def in_group(user, groups):
    """Returns a boolean if the user is in the given group, or comma-separated
    list of groups.

    Usage::

    {% if user|in_group:"Friends" %}
    ...
    {% endif %}

    or::

    {% if user|in_group:"Friends,Enemies" %}
    ...
    {% endif %}

    """
    group_list = groups.split(',')
    return bool(user.groups.filter(name__in=group_list).values('name'))
