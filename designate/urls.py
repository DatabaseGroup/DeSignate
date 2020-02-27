"""designate URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from webapp import views
from django.conf import settings
from django.conf.urls.static import static

ROOT_URLCONF = __name__

# general pattern:
# path(x, y, z) where ...
# x is the string thats added to the url, e.g if we call tree at 127.0.0.1:8000 x will add 'tree/' so the url is 127.0.0.1:8000/tree/ now
# y is the function that is called in the views.py file
# z is the name that is used to refer to this site (see base.html at e.g. {% url home %} for home)
urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.home), # adds home.html as goto file if the user doesnt pick home, tree, analysis or github
    path('home/', views.home, name="home"),
    path('files/', views.file_input, name="fileInput"),
    path('tree/', views.tree, name="tree"),
    path('names/', views.names, name="names"),
    path('results', views.results, name="results"),
    path('contact', views.contact, name="contact"),
    path('privacy', views.privacy, name="privacy"),
    path('impressum', views.impressum, name="impressum"),
    path('download_zip', views.download_zip, name="download_zip")
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
