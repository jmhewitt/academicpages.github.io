---
layout: archive
title: "Research"
permalink: /research/
author_profile: true
---

Please explore brief demonstrations of my research, mainly through (generally) 
reproducible code.  Refer to my [publications list](publications) or [CV](cv) for complete 
references and details about my research. 

{% for item in site.research %}
<h2>{{ item.title }}</h2>
<p><a href="{{ item.url }}">{{ item.description }}</a></p>
{% endfor %}
