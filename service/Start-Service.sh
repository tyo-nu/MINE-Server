source /vol/kbase/deployment/user-env.sh
source /vol/model-prod/anaconda3/bin/activate my-rdkit-env
uwsgi --master --processes 20 --cheaper 4 --http :7074 --http-timeout 600 --pidfile /vol/model-prod/mine-server/service/pid --daemonize /vol/model-prod/mine-server/service/error.log --wsgi-file /vol/model-prod/mine-server/lib/biokbase/mine_database/Server.py
