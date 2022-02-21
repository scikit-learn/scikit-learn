# Adapated from https://github.com/mxschmitt/action-tmate/blob/master/src/index.js
set -e

echo "Downloading tmate"
cd /tmp
wget https://github.com/tmate-io/tmate/releases/download/2.4.0/tmate-2.4.0-static-linux-amd64.tar.xz
tar xvf tmate-2.4.0-static-linux-amd64.tar.xz
tmate_path="/tmp/tmate-2.4.0-static-linux-amd64/tmate"
cd -

echo "Creating new session"

echo 'set +e' > /tmp/tmate.bashrc
set_default_command=(set-option -g default-command 'bash --rcfile /tmp/tmate.bashrc' ";")
# echo "number of opts = ${#set_default_command[@]}; opts are: ${set_default_command[@]}"

tmate_socket="/tmp/tmate.sock"
tmate="${tmate_path} -S ${tmate_socket}"

$tmate "${set_default_command[@]}" new-session -d
$tmate wait tmate-ready

echo "Created new session successfully"

echo "Fetching connection strings"
tmate_ssh=$(${tmate} display -p '#{tmate_ssh}')
tmate_web=$(${tmate} display -p '#{tmate_web}')

echo "Entering main loop"

while (( 1 )); do
    echo "Web shell: ${tmate_web}"
    echo "SSH: ${tmate_ssh}"

    sleep 5

    if [[ -e continue ]]; then
        echo "Found continue file, continuing"
        break
    fi

    if [[ ! -e "${tmate_socket}" ]]; then
        echo "tmate socket not found, exiting"
        break
    fi

done
