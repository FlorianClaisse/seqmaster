# Contig

## Find All

```bash
./Contig --findAll --inputA <path> --inputB <path> --type <nucl/prot > [--output <path>] [--accept <percentage>]
```

Permet à partir d'un fichier d'entrée au format fasta de déterminer qu'elles
Contig sont présent dans chaque fichier du dossier B.

Si le contig est présent
alors dans un fichier de sortie `output.txt` se trouvera le nom du fichier
concerné suis des contigs du fichier A trouvé. Puis un fichier
`<filename>-result.fasta` va se trouver le nom des contigs trouvé ainsi que la
sequance qui leur correspond.

De plus il est possible de définir un pourcentage `de 0 à 100` pour determiner
le pourcentage minimum de correspondance souhaité.

IMPORTANT ne pas lancer le programme dans un dossier qui contient deja des fichier de resultat.

## Codon Count

```bash
./Contig --codonCount --inputA <path> --output <path>
```
