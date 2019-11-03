#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals

AUTHOR = 'Oromion'
SITENAME = 'Finite element method'
SITEURL = 'http://localhost:8000'

PATH = 'content'

TIMEZONE = 'America/Lima'

DEFAULT_LANG = 'en'

ARTICLE_URL = 'blog/{date:%Y}/{date:%m}/{date:%d}/{slug}/'
ARTICLE_SAVE_AS = 'blog/{date:%Y}/{date:%m}/{date:%d}/{slug}/index.html'

THEME = 'themes/elegant'

MARKUP = ['md']
PLUGIN_PATHS = ['./plugins', './plugins_community']
PLUGINS = [
	'neighbors',
	'latex-prerender',
	'liquid_tags.img',
	'liquid_tags.video',
	'liquid_tags.include_code',
	'liquid_tags.literal',
	'pelican-ipynb.liquid',
	'pelican-ipynb.markup',
	'summary'
]

IGNORE_FILES = ['.ipynb_checkpoints']

# for liquid tags
CODE_DIR = './downloads/code'
NOTEBOOK_DIR = './downloads/notebooks'

IPYNB_USE_METACELL = True

STATIC_PATHS = ['downloads']

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
LINKS = (('Pelican', 'http://getpelican.com/'),
		 ('Python.org', 'http://python.org/'),
		 ('Jinja2', 'http://jinja.pocoo.org/'),
		 ('You can modify those links in your config file', '#'),)

# Social widget
SOCIAL = (('You can add links in your config file', '#'),
		  ('Another social link', '#'),)

DEFAULT_PAGINATION = 10

# Uncomment following line if you want document-relative URLs when developing
RELATIVE_URLS = True