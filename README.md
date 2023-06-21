# seqmaster

seqmaster is a program that contains several subroutines.
The main goal is to do different manipulations on genomic files to extract various relevant information.

Currently seqmaster has 3 sub-programs:
- findall : Lets you know which contig of an input file A is present in which file of a folder B.
- codonCount : Allows to calculate the codon bias.
- contigdiff :
- contigOrdered :


## Requirements

|         Requirement          | Version  |
|:----------------------------:|:--------:|
|  [GCC](https://gcc.gnu.org)  | >= 11.3  |
|  [CMake](https://cmake.org)  |  >= 3.4  |

## Dependencies

|                           Libs                            | Version  |
|:---------------------------------------------------------:|:--------:|
|   [seqan3](https://github.com/seqan/seqan3/tree/master)   | >= 1.1.0 |
|   [sharg-parser](https://github.com/seqan/sharg-parser)   | >= 3.2.0 |
| [csv-parser](https://github.com/vincentlaucsb/csv-parser) | >= 2.1.3 |

## Setup

```bash
git clone --recursive git@github.com:FlorianClaisse/seqmaster.git
cd seqmaster
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../sources
make
```

## Find All

```bash
./seqmaster findall -A <path> -tB <path> -t <nucl/prot> -o <path> [--accept <percentage>] [--threads <number>]
```

À partir d'un fichier d'entrée au format fasta détermine qu'elles Contig sont présent dans chaque fichier du dossier B.

Au niveau des paramètres de la ligne de commande, il y a :

- `-A --inputA` Le chemin vers le **fichier** (format fasta) contenant les contigs à trouver.
- `-B --inputB` Le chemin vers le **dossier** ou se trouve les fichiers à tester.
- `-t --type` Le type d'élément contenue dans les fichiers (nucl ou prot).
- `-o --output` Le chemin vers le **dossier** qui va contenir les fichiers de sortie.
- `--accept` Le pourcentage d'acceptation pour qu'un élément soit considéré comme reconnu par defaut il vaut 100%.
- `--threads` Le nombre de threads que le programme peut utiliser par défaut il vaut 4.

À partir du moment ou un Contig qui se trouve dans le fichier A est reconnu dans un fichier du dossier B,
alors dans le fichier `output.txt` va se trouver à la suite du nom du fichier de B en question le nom du contig
de A trouvé suivi du pourcentage si `--accept` a était donné dans la ligne de commande.
Le fichier `output.txt` est au format tabulé.

En plus du fichier `output.txt` pour chaque fichier du dossier d'entrée B un fichier `nomfichier-result.fasta`
est généré. Ce fichier se trouve dans le dossier d'output et contient les contigs trouvé dans le fichier de B
correspondant au format `nom-contig-B -> nom-contig-A -> pourcentage`. Le pourcentage n'est présent que dans le
cas ou il est spécifié en option de la ligne de commande. Puis un retour ligne est effectué afin de placer le bout
de contig reconnu. Dans le cas ou le pourcentage est different de 100, si `--type` vaut `nucl` c'est le bout du contig
présent dans B qui est donné sinon si `--type` vaut `prot` alors c'est toute la proteine de B qui est donné.

### Exemple d'output

- `output.txt`

```text
Filename
file1.fa    nameA -> percentage     name2A -> percentage
file2.fa
file3.fa    name23A -> percentage
```

- `file1-result.fasta`

```text
>Contig_1_B -> nameA -> percentage
contigValue
>Contig_5_B -> name2A -> percentage
```

## Codon Count

```bash
./seqmaster codonCount -A <path> -o <path>
```

À partir d'un fichier d'entrée, le programme va compter le nombre de chaque codon présent dans chaque
Contig du fichier d'entrée.

Le résultat de ce calcule se trouvera dans un fichier `output.txt` au format tabulé

### Exemple d'output

- `output.txt`

```text
Contig Name     Codon   Number  Percentage
>Contig_1       ATG     3       25%
>Contig_1       CTG     6       50%
>Contig_1       ATT     3       25%
>Contig_2       AAA     12      50%
>Contig_2       TGT     12      50%
```

## Contig Diff

```bash
./seqmaster contigdiff -A <path> -B <path> -t <nucl/prot> -o <path> --accept <percentage> --threads <num>
```

## Gene Mut

```bash
./seqmaster genemut -i <path> -g <path> -o <path>
```

```bash
./seqmaster contigordered -i <path> -o <path>
```
