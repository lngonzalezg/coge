#!/bin/bash
#------------------------------------------------------------------------------
# Script to start/stop Mojolicious development web server (morbo) for this sandbox.
#------------------------------------------------------------------------------

start() {
    export COGE_HOME=$(pwd)
    export IRODSENV=$COGE_HOME/irodsEnv_local
#    hypnotoad -f ./web/services/api.pl >> ./api.log 2>&1 &
#    ./web/services/api.pl daemon -l http://localhost:3304 >> ./api.log 2>&1 &
    port=$(grep MOJOLICIOUS_PORT ./coge.conf | cut -d ' ' -f 2)
    if [ "$port" == "" ]; then
        port=3303
    fi
    morbo -l http://localhost:$port ./web/services/api.pl
    echo "Started API (port $port)"
}

stop() {
    me=$(whoami)
    stopped=""

    # hypnotoad (prod mode): graceful stop via its pid file
    if [ -f /tmp/coge-api.pid ]; then
        hypnotoad ./web/services/api.pl --stop 2>/dev/null && stopped=1
    fi

    # morbo (dev mode): pgrep 'morbo' only matches the manager — the worker's
    # process name is api.pl, so it used to survive a stop and hold the port.
    # Kill both explicitly.
    pid=$(pgrep -u $me 'morbo')
    if [ "$pid" != "" ]; then
        kill -9 $pid
        stopped=1
    fi
    pids=$(pgrep -u $me -f 'web/services/api\.pl')
    if [ "$pids" != "" ]; then
        kill -9 $pids 2>/dev/null
        stopped=1
    fi

    if [ "$stopped" != "" ]; then
        echo "Stopped API"
    else
        echo "Not running"
    fi
}

prod() {
    # Production mode: hypnotoad preforks warm workers (no first-request module
    # load) and does not watch files (no accidental hot reloads). Runs in the
    # foreground for supervisor.
    export COGE_HOME=$(pwd)
    export IRODSENV=$COGE_HOME/irodsEnv_local
    exec hypnotoad -f ./web/services/api.pl
}

case "$1" in
    start)
        start
        exit 0
        ;;
    prod)
        prod
        ;;
    stop)
        stop
        exit 0
        ;;
    restart)
        stop
        sleep 1
        start
        exit 0
        ;;
    *)
        echo "Usage: $0 {start|prod|stop|restart}" 1>&2
        exit 1
        ;;
esac
