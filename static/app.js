async function* walkDir(dirHandle) {
  for await (const entry of dirHandle.values()) {
    if (entry.kind === 'file') {
      yield entry;
    } else if (entry.kind === 'directory') {
      yield* walkDir(entry);
    }
  }
}

async function pickAndUpload() {
  if (!window.showDirectoryPicker) {
    alert("Dein Browser unterstützt die File System Access API nicht. Bitte Chrome/Edge nutzen.");
    return;
  }
  const dir = await window.showDirectoryPicker();
  let files = [];
  for await (const handle of walkDir(dir)) {
    if (/\.(h5|csv)$/i.test(handle.name)) {
      files.push(handle);
    }
  }
  if (!files.length) {
    alert("Keine .h5 oder .csv gefunden.");
    return;
  }

  const prog = document.getElementById('prog');
  const status = document.getElementById('status');
  const formStart = document.getElementById('startLocal');

  prog.style.display = 'block';
  prog.max = files.length; prog.value = 0;
  status.textContent = `0 / ${files.length} Dateien`;

  // Parallele Uploads (limit 3)
  const limit = 3;
  let idx = 0, done = 0, errs = 0;

  async function uploadOne(handle) {
    const file = await handle.getFile();
    const form = new FormData();
    form.append('file', file, handle.name);
    const res = await fetch('/upload_fs', { method: 'POST', body: form });
    if (!res.ok) throw new Error(await res.text());
  }

  async function worker() {
    while (true) {
      let myIdx;
      if (idx < files.length) {
        myIdx = idx++;
      } else break;
      try {
        await uploadOne(files[myIdx]);
      } catch (e) {
        console.error(e); errs++;
      }
      done++;
      prog.value = done;
      status.textContent = `${done} / ${files.length} Dateien hochgeladen${errs? " ("+errs+" Fehler)": ""}`;
    }
  }

  const workers = Array.from({length: limit}, worker);
  await Promise.all(workers);

  if (errs) {
    status.innerHTML += `<div class="text-danger mt-2">${errs} Upload-Fehler</div>`;
  } else {
    status.innerHTML += `<div class="text-success mt-2">Upload abgeschlossen.</div>`;
  }
  formStart.style.display = 'block';
}

document.addEventListener('DOMContentLoaded', () => {
  const btn = document.getElementById('pick');
  if (btn) btn.addEventListener('click', pickAndUpload);
});
