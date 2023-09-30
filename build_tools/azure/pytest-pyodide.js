const { opendir } = require('node:fs/promises');
const { loadPyodide } = require("pyodide");

async function main() {
    let exit_code = 0;
    try {
        global.pyodide = await loadPyodide();
        let pyodide = global.pyodide;
        const FS = pyodide.FS;
        const NODEFS = FS.filesystems.NODEFS;

        let mountDir = "/mnt";
        pyodide.FS.mkdir(mountDir);
        pyodide.FS.mount(pyodide.FS.filesystems.NODEFS, { root: "." }, mountDir);

        await pyodide.loadPackage(["micropip"]);
        await pyodide.runPythonAsync(`
            import glob
            import micropip

            wheels = glob.glob('/mnt/dist/*.whl')
            wheels = [f'emfs://{wheel}' for wheel in wheels]
            print(f'installing wheels: {wheels}')
            await micropip.install(wheels);

            pkg_list = micropip.list()
            print(pkg_list)
        `);

        // Pyodide is built without OpenMP, need to set environment variable to
        // skip related test
        await pyodide.runPythonAsync(`
            import os
            os.environ['SKLEARN_SKIP_OPENMP_TEST'] = 'true'
        `);

        await pyodide.runPythonAsync("import micropip; micropip.install('pytest')");
        let pytest = pyodide.pyimport("pytest");
        let args = process.argv.slice(2);
        console.log('pytest args:', args);
        exit_code = pytest.main(pyodide.toPy(args));
    } catch (e) {
        console.error(e);
        // Arbitrary exit code here. I have seen this code reached instead of a
        // Pyodide fatal error sometimes
        exit_code = 66;

    } finally {
        process.exit(exit_code);
    }
}

main();
