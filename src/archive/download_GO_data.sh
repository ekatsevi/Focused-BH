# Download GO ontology and gene annotations from the web

data_dir="/home/ekatsevi/Local/large_file_storage/filtered_BH/DAG/data"
ontology_file=$data_dir"/go-basic.obo"
annotations_zipped_file=$data_dir"/goa_human.gaf.gz"
annotations_file=$data_dir"/goa_human.gaf"

if test -f $ontology_file; then
    echo "Ontology file already downloaded!"
else 
    echo "Downloading ontology file..."
    wget http://current.geneontology.org/ontology/go-basic.obo -P $data_dir    
fi

if test -f $annotations_file; then
    echo "Annotations file already downloaded and unzipped!"
else
    if test -f $annotations_zipped_file; then
    echo "Annotations file already downloaded!"
    echo "Unzipping annotations file..."
    gunzip $annotations_zipped_file
    else
        echo "Downloading annotations file..."
        wget http://current.geneontology.org/annotations/goa_human.gaf.gz -P $data_dir
        echo "Unzipping annotations file..."
        gunzip $annotations_zipped_file
    fi
fi
echo "Download complete."