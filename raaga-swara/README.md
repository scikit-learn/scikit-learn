# Raaga Swara

A small, dependency-free browser app for recording raaga examples and practicing
recognition by listening.

## Run it

Microphone access requires a secure browser context. From the repository root,
start a local web server:

```sh
python -m http.server 8000
```

Then open <http://localhost:8000/raaga-swara/> in a browser that supports
`MediaRecorder`, `getUserMedia`, and IndexedDB. Allow microphone access when
prompted, enter a raaga name, and start recording. Each recording stops after
30 seconds and is saved to the local repository automatically.

## Features

- Record a 30-second clip from your microphone and label it with a
  raaga name. Suggested names are available, and custom names are accepted.
- Listen to or delete saved recordings in the repository.
- Start a listening quiz after saving recordings for at least two different
  raagas. Answer choices are drawn from the names in your repository.
- Clips persist in this browser's IndexedDB database and are not uploaded.

The repository is local to the browser profile and device. Clearing browser site
data removes the recordings; this app does not currently sync or back up clips
between devices.
