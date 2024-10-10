export NXF_VER=22.10.8

if [ -f "/usr/local/bin/nextflow" ]; then
    echo "Nextflow is already installed."
else
    wget https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow
    bash nextflow
    chmod +x nextflow
    mv nextflow /usr/local/bin
    echo "Nextflow v$NXF_VER installed successfully."
fi