psql --host=ccpglsstdev.in2p3.fr --port=6553 --username=$(whoami) --dbname=$(whoami)

GRANT ALL PRIVILEGES ON SCHEMA main_20230111 TO boutigny;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA main_20230111 TO boutigny;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA main_20230111 TO boutigny;
