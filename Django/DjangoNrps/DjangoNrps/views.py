from django.http import HttpResponse

def index(request):
    return HttpResponse("Hello, you're at the iGEM Team Heidelberg's NRPS Designer index.")

