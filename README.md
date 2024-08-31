# code

###############################################################################################
##### Les différentes modifications acec le modèle avec terme (x-a_c)^+(I_t-mean(I_t))^+kt5
###############################################################################################

### Intervention sur les méthodes de simulation des intervalles de confiance du modèle.

# En utilisant la méthode dans PLAT avec des ordres de ARIMA appropriés (puisque la remarque avec PLAT c'est que les choix d'ordre Arima n'était pas basé sur la fonction auto.arima mais plutôt avec des à priori permettant d'éviter des problèmes dans les prévisions), j'arrive à reproduire approximativement les résultats de PLAT, mais cette méthode appliquée à notre modèle fait n'importe quoi dans les prévisions (en gros les prévisions sont très volatiles tout âge confondu.)

### En reprenant la méthode boostrap (avec la technique utilisée dans PLAT c'est-à-dire le choix de 2,5% des données en sus et en sous), on tombe sur des résultats un peu mieux avec notre modèle c'est-à-dire moins volatile et ceci, seulement si l'ordre de simulation est dans les 100 ou 200, après cela ça fait du n'importe quoi.

### Avec un auto.arima sur les kt on arrive à quelque chose mais à partir d'un certain temps, le taux de mortalité se redresse sur les intervalles de confiance.

### Par contre lorsqu'on essaie PLAT avec gamma simulé normalement sur toute la période de cohorte, on a de meilleur résultat pour l'IC.


### J'ai aussi essayé de voir le cas où je fite les paramètres temporelles avec des Arima comme dans PLAT, mais là également on remarque une divergence de l'IC

### En essayant la simulation sur tous les autres paramètres du modèle et pas sur kt5, on observe toujours un problème de redressement dans l'IC.


###########################################################################
### En utilisant le boostrap pour B petit (100) et en utilisant les quantiles
### d'ordre 1-a/2 pour calculer l'IC, on se retrouve finalement avec IC qui contient 
### nos prévisions



###### Grosse remarque
En réexaminant les simulations des trajectoires futures obtenues avec le package StMoMo pour le modèle utilisé dans le mémoire de Pape Falh, on se rend compte que la simulation pour obtenir l'intervalle de confiance part de l'hypothèse que les "Kappa_t" suivent une marche aléatoire multivariée, en changeant cette option en Arima indépendant, on se rend compte que l'intervalle de confiance obtenue diverge, exactement comme dans mon cas.






