env_push() {
    eval value=\$$1
    if [[ $value = "" ]]; then
        export $1=$2
    elif [ -d "$2" ]; then
        tmp=$(echo $value | tr ':' '\n' | awk '$0 != "'"$2"'"' | paste -sd: -)
        if [[ $tmp = "" ]]; then export $1=$2; else export $1=$2:$tmp; fi
    fi
}
source /vol/kbase/deployment/user-env.sh
env_push PATH /vol/model-prod/Python2.7.2/bin
env_push PYTHONPATH /vol/model-prod/mine-server/lib
env_push PYTHONPATH /vol/model-prod/Python2.7.2/lib
env_push PYTHONPATH /vol/model-prod/Python2.7.2/lib/python2.7
env_push PYTHONPATH /vol/model-prod/Python2.7.2/lib/python2.7/site-packages
env_push LD_LIBRARY_PATH /vol/model-prod/Python2.7.2/lib/python2.7/site-packages/openbabel/lib
export PATH
export PYTHONPATH
echo $PYTHONPATH
which python
export LD_LIBRARY_PATH
uwsgi --master --processes 20 --cheaper 4 --http :7074 --http-timeout 600 --pidfile /vol/model-prod/mine-server/service/pid --daemonize /vol/model-prod/mine-server/service/error.log --wsgi-file /vol/model-prod/mine-server/lib/biokbase/mine_database/Server.py
