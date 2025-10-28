/* define several functions to replace jQuery methods
 * inspired by https://tobiasahlin.com/blog/move-from-jquery-to-vanilla-javascript/
 */

/**
 * Execute a method if DOM has finished loading
 *
 * @param {function} callback the method to execute
 */
export function documentReady(callback) {
  if (document.readyState != "loading") callback();
  else document.addEventListener("DOMContentLoaded", callback);
}
