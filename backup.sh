rsync -avh --progress --partial --delete \
  --exclude='.git/' \
  --exclude='.venv/' \
  --exclude='__pycache__/' \
  --exclude='.Rproj.user/' \
  --exclude='.Rhistory' \
  --exclude='.RData' \
  --exclude='.DS_Store' \
  --exclude='Thumbs.db' \
  ~/work/qPCR/ \
  /Volumes/nipper/qPCR/

