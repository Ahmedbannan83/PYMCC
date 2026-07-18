# Calcul de TPR20,10 et kQ pour l'étalonnage absolu CyberKnife

## Objectif

Le module `pymcc.cyberknife_calibration` aide à déterminer :

1. le `TPR20,10` équivalent d'un faisceau CyberKnife mesuré avec le cône fixe de 60 mm ;
2. le facteur de qualité `kQ` correspondant à une chambre d'ionisation sélectionnée.

Il applique la procédure de conversion décrite dans le chapitre 2 du *Physics Essentials Guide* Accuray, avec les données BJR Supplement 25 et la Table 16 du TRS-398 (Rev. 1).

> **Avertissement clinique**  
> Cet outil est une aide de calcul destinée à un physicien médical qualifié. Il ne remplace pas un protocole d'étalonnage, une revue indépendante, ni les exigences locales de QA.

## Géométrie requise par la procédure Accuray

La procédure d'étalonnage absolu demande une mesure PDD avec :

- collimateur fixe : **60 mm** ;
- SSD : **1000 mm** ;
- mesure à 100 mm et 200 mm de profondeur ;
- normalisation à (d_{max}).

Le cône circulaire de 60 mm, défini à 800 mm SAD, est converti en carré équivalent à SSD = 1000 mm :

[
ESFS_{1000} = 1.125 	imes 60 	ext{mm} = 67.5 	ext{mm}.
]

Le fichier `PDD.mcc` actuellement présent dans le dépôt indique `SSD=800 mm`. Le module le détecte et retourne un avertissement. Une valeur obtenue avec cette géométrie ne doit pas être utilisée pour l'étalonnage clinique sans mesure conforme à SSD = 1000 mm.

## Étapes du calcul

1. `calculate_from_mcc()` lit le premier PDD photonique du fichier MCC et extrait :
   - `PDD100` : dose relative à 100 mm ;
   - `PDD200` : dose relative à 200 mm.

2. Les valeurs BJR-25, à 6 cm, 7 cm et 10 cm de champ carré, sont interpolées linéairement en taille de champ à **6.75 cm**.

3. Le `PDD100` mesuré est comparé aux valeurs BJR-25 afin d'inférer l'énergie équivalente dans la plage de référence 4-8 MV.

4. La correction de taille de champ BJR-25 est appliquée séparément aux PDD à 100 mm et 200 mm afin d'obtenir les PDD équivalents à **10 cm x 10 cm**.

5. Le module emploie la relation déjà utilisée dans `pymcc.wtscans.PDD.calc_q_index()` :

[
TPR_{20,10} = 1.2661 	imes rac{PDD_{200}^{eq}}{PDD_{100}^{eq}} - 0.0595.
]

6. Si un modèle de chambre est fourni, le code cherche le modèle dans la Table 16 du TRS-398 et interpole linéairement `kQ` en fonction du `TPR20,10` calculé.

## Exemple

```python
from pymcc.cyberknife_calibration import calculate_from_mcc, result_as_dict

result = calculate_from_mcc(
    "PDD.mcc",
    ionization_chamber="PTW 31021 Semiflex 3D",
)

print(result_as_dict(result))
```

Le résultat contient notamment :

- `equivalent_square_mm_at_1000_ssd` : 67.5 mm pour le cône 60 mm ;
- `bjr25_inferred_energy_mv` ;
- `equivalent_pdd100_10x10` et `equivalent_pdd200_10x10` ;
- `tpr20_10` ;
- `ionization_chamber` et `kq_trs398_table16` ;
- `warnings`, à examiner avant tout usage.

## Modèles de chambres

Les modèles de la Table 16 incluent notamment :

- PTW 30010, 30011, 30012, 30013 ;
- PTW 31010, 31013, 31016, 31021 et 31022 ;
- Exradin A1SL, A12, A12S, A18, A19, A26 et A28 ;
- IBA CC13, CC25, FC23-C, FC65-G et FC65-P ;
- NE 2561/2611A et NE 2571 ;
- Capintec PR-06C ;
- Sun Nuclear SNC125c et SNC600c.

Les données de la Table 16 couvrent `TPR20,10 = 0.56` à `0.82`. Le code refuse une extrapolation hors de cette plage.

## Références

1. Accuray Incorporated. *Physics Essentials Guide*, 1075879-ENG A, Chapter 2, **Absolute LINAC Dose Calibration**, pp. 2-155 à 2-159.  
   Procédure pour la mesure avec le cône fixe 60 mm, conversion en carré équivalent 67.5 mm et comparaison avec BJR.

2. Jordan TJ. *Megavoltage X-ray Beams: 2-50 MV*. British Journal of Radiology Supplement 25, 1996.  
   Tables de PDD à SSD = 100 cm utilisées pour l'interpolation de taille de champ et l'équivalence 10 cm x 10 cm.

3. International Atomic Energy Agency. *Absorbed Dose Determination in External Beam Radiotherapy: An International Code of Practice for Dosimetry Based on Standards of Absorbed Dose to Water*. IAEA TRS-398 (Rev. 1), 2024, **Table 16**, pp. 90-93.  
   Valeurs de `kQ` pour chambres cylindriques en fonction de `TPR20,10`.

4. Almond PR, Biggs PJ, Coursey BM, et al. AAPM's TG-51 protocol for clinical reference dosimetry of high-energy photon and electron beams. *Medical Physics*. 1999;26(9):1847-1870.  
   Relation de qualité de faisceau utilisée par le module.

## Limites

- La conversion BJR-25 utilise les tables de référence 4, 5, 6 et 8 MV ; le code n'extrapole pas au-delà de cette plage.
- Les mesures PDD doivent être revues pour la géométrie, la profondeur, le détecteur, la stabilité et la correction de volume avant emploi clinique.
- Le module ne calcule pas les autres corrections de l'étalonnage absolu : (k_{TP}), (k_{pol}), (k_s), correction d'électromètre, (N_{D,w}), ni le facteur de volume moyen `Prp` éventuel pour les chambres Farmer en faisceau FFF.
