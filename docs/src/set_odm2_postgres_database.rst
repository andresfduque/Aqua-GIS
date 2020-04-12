1. sudo -i -s postgres
2. psql
3. CREATE DATABASE odm2col;
4. \q
5. psql -d odm2col -a -f odm2SQLFile.sql
6. python CvLoad.py postgresql+psycopg2://postgres:postgres@localhost/odm2col [username:password@host/dbname]
6. python Odm2ColStructure.py postgresql+psycopg2://postgres:postgres@localhost/odm2col [username:password@host/dbname]
