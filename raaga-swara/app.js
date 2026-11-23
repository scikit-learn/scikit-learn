"use strict";

const RECORDING_DURATION_SECONDS = 30;
const DATABASE_NAME = "raaga-swara-library";
const DATABASE_VERSION = 1;
const CLIP_STORE = "clips";

const recordButton = document.querySelector("#record-button");
const recordStatus = document.querySelector("#record-status");
const recordProgress = document.querySelector("#record-progress");
const ragaInput = document.querySelector("#raga-input");
const quizCount = document.querySelector("#quiz-count");
const quizContent = document.querySelector("#quiz-content");
const libraryCount = document.querySelector("#library-count");
const libraryList = document.querySelector("#library-list");

let database;
let clips = [];
let mediaRecorder;
let recordingStream;
let recordingTimer;
let recordingStartedAt;
let recordingRaga;
let lastCountdownSecond;
let recordingFailed = false;
let recordingChunks = [];
let libraryObjectUrls = [];
let quizObjectUrl;
let currentQuizClip;
let quizDeck = [];

function openDatabase() {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(DATABASE_NAME, DATABASE_VERSION);

    request.onupgradeneeded = () => {
      const store = request.result.createObjectStore(CLIP_STORE, {
        keyPath: "id",
      });
      store.createIndex("createdAt", "createdAt");
      store.createIndex("raga", "raga");
    };
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
    request.onblocked = () => reject(new Error("The raaga library is busy in another tab. Close that tab and try again."));
  });
}

function saveClip(clip) {
  return new Promise((resolve, reject) => {
    const transaction = database.transaction(CLIP_STORE, "readwrite");
    transaction.objectStore(CLIP_STORE).add(clip);
    transaction.oncomplete = () => resolve();
    transaction.onerror = () => reject(transaction.error);
    transaction.onabort = () => reject(transaction.error || new Error("The recording could not be saved."));
  });
}

function readClips() {
  return new Promise((resolve, reject) => {
    const transaction = database.transaction(CLIP_STORE, "readonly");
    const request = transaction.objectStore(CLIP_STORE).getAll();
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
}

function removeClip(id) {
  return new Promise((resolve, reject) => {
    const transaction = database.transaction(CLIP_STORE, "readwrite");
    transaction.objectStore(CLIP_STORE).delete(id);
    transaction.oncomplete = () => resolve();
    transaction.onerror = () => reject(transaction.error);
    transaction.onabort = () => reject(transaction.error || new Error("The clip could not be deleted."));
  });
}

function setRecordStatus(message, state = "") {
  recordStatus.textContent = message;
  recordStatus.classList.toggle("is-error", state === "error");
  recordStatus.classList.toggle("is-success", state === "success");
}

function getRecordingMimeType() {
  if (!MediaRecorder.isTypeSupported) {
    return "";
  }

  return [
    "audio/webm;codecs=opus",
    "audio/webm",
    "audio/mp4",
  ].find((mimeType) => MediaRecorder.isTypeSupported(mimeType)) || "";
}

function makeClipName(raga) {
  return `${raga} · ${new Date().toLocaleString()}`;
}

async function startRecording() {
  const raga = ragaInput.value.trim();
  if (!raga) {
    ragaInput.focus();
    setRecordStatus("Enter a raaga name before recording.", "error");
    return;
  }
  if (!navigator.mediaDevices?.getUserMedia || !window.MediaRecorder) {
    setRecordStatus("Audio recording is not supported in this browser.", "error");
    return;
  }

  recordButton.disabled = true;
  setRecordStatus("Requesting microphone access…");
  try {
    recordingStream = await navigator.mediaDevices.getUserMedia({ audio: true });
    recordingRaga = raga;
    recordingChunks = [];
    recordingFailed = false;
    const mimeType = getRecordingMimeType();
    mediaRecorder = new MediaRecorder(
      recordingStream,
      mimeType ? { mimeType } : undefined,
    );
    mediaRecorder.addEventListener("dataavailable", (event) => {
      if (event.data.size > 0) {
        recordingChunks.push(event.data);
      }
    });
    mediaRecorder.addEventListener("stop", finishRecording, { once: true });
    mediaRecorder.addEventListener("error", (event) => {
      const message = event.error?.message || "Unknown media recorder error.";
      recordingFailed = true;
      setRecordStatus(`Recording failed: ${message}`, "error");
    }, { once: true });
    mediaRecorder.start();
    recordingStartedAt = Date.now();
    lastCountdownSecond = RECORDING_DURATION_SECONDS;
    recordButton.classList.add("is-recording");
    recordButton.innerHTML = '<span class="record-dot" aria-hidden="true"></span> Recording…';
    setRecordStatus(`Recording for ${recordingRaga}. It will stop automatically at 30 seconds.`);
    recordProgress.value = 0;
    recordingTimer = window.setInterval(updateRecordingProgress, 100);
    window.setTimeout(() => {
      if (mediaRecorder?.state === "recording") {
        mediaRecorder.stop();
      }
    }, RECORDING_DURATION_SECONDS * 1000);
  } catch (error) {
    stopRecordingStream();
    recordButton.disabled = false;
    setRecordStatus(`Could not start recording: ${error.message}`, "error");
  }
}

function updateRecordingProgress() {
  const elapsedSeconds = Math.min(
    (Date.now() - recordingStartedAt) / 1000,
    RECORDING_DURATION_SECONDS,
  );
  recordProgress.value = elapsedSeconds;
  const remaining = Math.ceil(RECORDING_DURATION_SECONDS - elapsedSeconds);
  if (remaining !== lastCountdownSecond) {
    lastCountdownSecond = remaining;
    setRecordStatus(`Recording for ${recordingRaga}. ${remaining} seconds remaining.`);
  }
}

async function finishRecording() {
  window.clearInterval(recordingTimer);
  recordProgress.value = RECORDING_DURATION_SECONDS;
  recordButton.disabled = true;
  recordButton.classList.remove("is-recording");
  recordButton.innerHTML = '<span class="record-dot" aria-hidden="true"></span> Start recording';
  stopRecordingStream();

  if (recordingFailed) {
    recordingChunks = [];
    recordProgress.value = 0;
    recordButton.disabled = false;
    return;
  }

  const mimeType = mediaRecorder.mimeType || recordingChunks[0]?.type || "audio/webm";
  const blob = new Blob(recordingChunks, { type: mimeType });
  if (blob.size === 0) {
    setRecordStatus("The recording was empty and could not be saved. Please try again.", "error");
    recordProgress.value = 0;
    recordButton.disabled = false;
    return;
  }

  const clip = {
    id: crypto.randomUUID(),
    raga: recordingRaga,
    name: makeClipName(recordingRaga),
    createdAt: new Date().toISOString(),
    durationMs: Date.now() - recordingStartedAt,
    mimeType,
    blob,
  };

  setRecordStatus("Saving your clip to this browser…");
  try {
    await saveClip(clip);
    await refreshLibrary();
    setRecordStatus(`Saved a ${formatDuration(clip.durationMs)} clip of ${clip.raga}.`, "success");
  } catch (error) {
    setRecordStatus(`Could not save your recording: ${error.message}`, "error");
  } finally {
    recordingChunks = [];
    recordProgress.value = 0;
    recordingFailed = false;
    recordButton.disabled = false;
  }
}

function stopRecordingStream() {
  if (recordingStream) {
    recordingStream.getTracks().forEach((track) => track.stop());
    recordingStream = undefined;
  }
}

function formatDuration(durationMs) {
  const seconds = Math.round(durationMs / 1000);
  return `${seconds} sec`;
}

function revokeLibraryObjectUrls() {
  libraryObjectUrls.forEach((url) => URL.revokeObjectURL(url));
  libraryObjectUrls = [];
}

function renderLibrary() {
  revokeLibraryObjectUrls();
  libraryCount.textContent = `${clips.length} ${clips.length === 1 ? "clip" : "clips"}`;
  if (clips.length === 0) {
    libraryList.innerHTML = '<p class="empty-state">No recordings yet. Your saved clips will appear here.</p>';
    return;
  }

  libraryList.replaceChildren();
  for (const clip of clips) {
    const card = document.createElement("article");
    card.className = "clip-card";
    const details = document.createElement("div");
    details.className = "clip-details";
    const title = document.createElement("p");
    title.className = "clip-title";
    title.textContent = clip.raga;
    const meta = document.createElement("p");
    meta.className = "clip-meta";
    meta.textContent = `${new Date(clip.createdAt).toLocaleString()} · ${formatDuration(clip.durationMs)}`;
    details.append(title, meta);

    const deleteButton = document.createElement("button");
    deleteButton.className = "delete-button";
    deleteButton.type = "button";
    deleteButton.textContent = "Delete";
    deleteButton.setAttribute("aria-label", `Delete ${clip.raga} recording`);
    deleteButton.addEventListener("click", () => deleteSavedClip(clip));

    const audio = document.createElement("audio");
    audio.className = "clip-audio";
    audio.controls = true;
    audio.preload = "none";
    const url = URL.createObjectURL(clip.blob);
    libraryObjectUrls.push(url);
    audio.src = url;
    card.append(details, deleteButton, audio);
    libraryList.append(card);
  }
}

async function deleteSavedClip(clip) {
  if (!window.confirm(`Delete this recording of ${clip.raga}?`)) {
    return;
  }
  try {
    await removeClip(clip.id);
    await refreshLibrary();
  } catch (error) {
    setRecordStatus(`Could not delete the clip: ${error.message}`, "error");
  }
}

async function refreshLibrary() {
  clips = (await readClips()).sort(
    (first, second) => second.createdAt.localeCompare(first.createdAt),
  );
  const savedClipIds = new Set(clips.map((clip) => clip.id));
  quizDeck = quizDeck.filter((clip) => savedClipIds.has(clip.id));
  renderLibrary();
  updateQuizAvailability();
}

function updateQuizAvailability() {
  const ragas = [...new Set(clips.map((clip) => clip.raga))].sort((a, b) => a.localeCompare(b));
  quizCount.textContent = ragas.length < 2
    ? `${ragas.length} of 2 raagas`
    : `${ragas.length} raagas`;

  if (ragas.length < 2) {
    quizContent.innerHTML = '<p class="empty-state">Your quiz will be ready when your library has clips from at least two raagas.</p>';
    currentQuizClip = undefined;
    return;
  }

  if (!currentQuizClip || !clips.some((clip) => clip.id === currentQuizClip.id)) {
    startQuizRound(ragas);
  } else {
    renderQuiz(ragas);
  }
}

function drawQuizClip() {
  if (quizDeck.length === 0) {
    quizDeck = [...clips];
    for (let index = quizDeck.length - 1; index > 0; index -= 1) {
      const swapIndex = Math.floor(Math.random() * (index + 1));
      [quizDeck[index], quizDeck[swapIndex]] = [quizDeck[swapIndex], quizDeck[index]];
    }
  }
  return quizDeck.pop();
}

function startQuizRound(ragas) {
  currentQuizClip = drawQuizClip();
  renderQuiz(ragas);
}

function renderQuiz(ragas) {
  if (quizObjectUrl) {
    URL.revokeObjectURL(quizObjectUrl);
  }
  quizObjectUrl = URL.createObjectURL(currentQuizClip.blob);
  quizContent.replaceChildren();

  const audio = document.createElement("audio");
  audio.className = "quiz-audio";
  audio.controls = true;
  audio.preload = "metadata";
  audio.src = quizObjectUrl;

  const form = document.createElement("form");
  form.className = "quiz-form";
  const label = document.createElement("label");
  label.htmlFor = "quiz-answer";
  label.textContent = "Which raaga do you hear?";
  const select = document.createElement("select");
  select.id = "quiz-answer";
  select.required = true;
  const initialOption = document.createElement("option");
  initialOption.value = "";
  initialOption.textContent = "Choose a raaga…";
  initialOption.disabled = true;
  initialOption.selected = true;
  select.append(initialOption);
  ragas.forEach((raga) => {
    const option = document.createElement("option");
    option.value = raga;
    option.textContent = raga;
    select.append(option);
  });
  const submitButton = document.createElement("button");
  submitButton.className = "button button-primary";
  submitButton.type = "submit";
  submitButton.textContent = "Check answer";
  form.append(label, select, submitButton);
  form.addEventListener("submit", (event) => {
    event.preventDefault();
    if (!select.value) {
      return;
    }
    const isCorrect = select.value === currentQuizClip.raga;
    feedback.textContent = isCorrect
      ? `Correct — this clip is ${currentQuizClip.raga}.`
      : `Not quite — this clip is ${currentQuizClip.raga}.`;
    feedback.classList.toggle("is-correct", isCorrect);
    feedback.classList.toggle("is-wrong", !isCorrect);
    submitButton.disabled = true;
    select.disabled = true;
    nextButton.hidden = false;
  });

  const feedback = document.createElement("p");
  feedback.className = "feedback";
  feedback.setAttribute("role", "status");
  const nextButton = document.createElement("button");
  nextButton.className = "button button-secondary";
  nextButton.type = "button";
  nextButton.textContent = "Next clip";
  nextButton.hidden = true;
  nextButton.addEventListener("click", () => startQuizRound(ragas));

  quizContent.append(audio, form, feedback, nextButton);
}

async function initialize() {
  try {
    database = await openDatabase();
    await refreshLibrary();
  } catch (error) {
    setRecordStatus(`Could not open your local raaga repository: ${error.message}`, "error");
    recordButton.disabled = true;
    libraryList.innerHTML = '<p class="empty-state">The local library is unavailable in this browser.</p>';
  }
}

recordButton.addEventListener("click", startRecording);
window.addEventListener("beforeunload", () => {
  window.clearInterval(recordingTimer);
  if (mediaRecorder?.state === "recording") {
    mediaRecorder.stop();
  }
  stopRecordingStream();
});

initialize();
