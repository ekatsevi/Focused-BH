# Download GO ontology and gene annotations from the web

base_dir=$1
data_dir=$base_dir"/data/raw/GO"
ontology_file=$data_dir"/go-basic.obo"
annotations_zipped_file=$data_dir"/goa_human.gaf.gz"
annotations_file=$data_dir"/goa_human.gaf"
gene_list_file=$data_dir"/gene_list.txt"

if test ! -f $ontology_file; then
    echo "Downloading ontology file..."
    wget http://current.geneontology.org/ontology/go-basic.obo -P $data_dir    
fi

if test ! -f $annotations_file; then
    if test -f $annotations_zipped_file; then
      gunzip $annotations_zipped_file
    else
        echo "Downloading annotations file..."
        wget http://current.geneontology.org/annotations/goa_human.gaf.gz -P $data_dir
        echo "Unzipping annotations file..."
        gunzip $annotations_zipped_file
    fi
fi

if test ! -f $gene_list_file; then
    echo "Downloading Van't'veer gene list..."
    wget http://cbl-gorilla.cs.technion.ac.il/VantVeerMoreLess5.txt -P $data_dir -O gene_list.txt
fi

# echo "Download complete."