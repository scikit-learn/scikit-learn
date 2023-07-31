"""
The :mod:`sklearn.utils.themes` module provides a
CSS template class for creating dynamic stylesheets.
It exposes utilities and themes for styling displayed diagrams.
"""

AUTO = """
#$id {
  --sklearn-color-1: black;
  --sklearn-color-2: #f0f8ff;
  --sklearn-color-3: #696969;
  --sklearn-color-4: white;
  --sklearn-color-5: #d4ebff;
  --sklearn-color-6: gray;

  @media (prefers-color-scheme: dark) {
    --sklearn-color-1: white;
    --sklearn-color-2: #424242;
    --sklearn-color-3: #696969;
    --sklearn-color-4: black;
    --sklearn-color-5: #616161;
    --sklearn-color-6: #616161;
  }
}

#$id {
  color: var(--sklearn-color-1);
}
#$id pre{
  padding: 0;
}
#$id div.sk-toggleable {
  background-color: var(--sklearn-color-4);
}
#$id label.sk-toggleable__label {
  cursor: pointer;
  display: block;
  width: 100%;
  margin-bottom: 0;
  padding: 0.3em;
  box-sizing: border-box;
  text-align: center;
}
#$id label.sk-toggleable__label-arrow:before {
  content: "▸";
  float: left;
  margin-right: 0.25em;
  color: var(--sklearn-color-3);
}
#$id label.sk-toggleable__label-arrow:hover:before {
  color: var(--sklearn-color-1);
}
#$id div.sk-estimator:hover label.sk-toggleable__label-arrow:before {
  color: var(--sklearn-color-1);
}
#$id div.sk-toggleable__content {
  max-height: 0;
  max-width: 0;
  overflow: hidden;
  text-align: left;
  background-color: var(--sklearn-color-2);
}
#$id div.sk-toggleable__content pre {
  margin: 0.2em;
  color: var(--sklearn-color-1);
  border-radius: 0.25em;
  background-color: var(--sklearn-color-2);
}
#$id input.sk-toggleable__control:checked~div.sk-toggleable__content {
  max-height: 200px;
  max-width: 100%;
  overflow: auto;
}
#$id input.sk-toggleable__control:checked~label.sk-toggleable__label-arrow:before {
  content: "▾";
}
#$id div.sk-estimator input.sk-toggleable__control:checked~label.sk-toggleable__label {
  background-color: var(--sklearn-color-5);
}
#$id div.sk-label input.sk-toggleable__control:checked~label.sk-toggleable__label {
  background-color: var(--sklearn-color-5);
}
#$id input.sk-hidden--visually {
  border: 0;
  clip: rect(1px 1px 1px 1px);
  clip: rect(1px, 1px, 1px, 1px);
  height: 1px;
  margin: -1px;
  overflow: hidden;
  padding: 0;
  position: absolute;
  width: 1px;
}
#$id div.sk-estimator {
  font-family: monospace;
  background-color: var(--sklearn-color-2);
  border: 1px dotted var(--sklearn-color-1);
  border-radius: 0.25em;
  box-sizing: border-box;
  margin-bottom: 0.5em;
}
#$id div.sk-estimator:hover {
  background-color: var(--sklearn-color-5);
}
#$id div.sk-parallel-item::after {
  content: "";
  width: 100%;
  border-bottom: 1px solid var(--sklearn-color-6);
  flex-grow: 1;
}
#$id div.sk-label:hover label.sk-toggleable__label {
  background-color: var(--sklearn-color-5);
}
#$id div.sk-serial::before {
  content: "";
  position: absolute;
  border-left: 1px solid var(--sklearn-color-6);
  box-sizing: border-box;
  top: 0;
  bottom: 0;
  left: 50%;
  z-index: 0;
}
#$id div.sk-serial {
  display: flex;
  flex-direction: column;
  align-items: center;
  background-color: var(--sklearn-color-4);
  padding-right: 0.2em;
  padding-left: 0.2em;
  position: relative;
}
#$id div.sk-item {
  position: relative;
  z-index: 1;
}
#$id div.sk-parallel {
  display: flex;
  align-items: stretch;
  justify-content: center;
  background-color: var(--sklearn-color-4);
  position: relative;
}
#$id div.sk-item::before, #$id div.sk-parallel-item::before {
  content: "";
  position: absolute;
  border-left: 1px solid var(--sklearn-color-6);
  box-sizing: border-box;
  top: 0;
  bottom: 0;
  left: 50%;
  z-index: -1;
}
#$id div.sk-parallel-item {
  display: flex;
  flex-direction: column;
  z-index: 1;
  position: relative;
  background-color: var(--sklearn-color-4);
}
#$id div.sk-parallel-item:first-child::after {
  align-self: flex-end;
  width: 50%;
}
#$id div.sk-parallel-item:last-child::after {
  align-self: flex-start;
  width: 50%;
}
#$id div.sk-parallel-item:only-child::after {
  width: 0;
}
#$id div.sk-dashed-wrapped {
  border: 1px dashed var(--sklearn-color-6);
  margin: 0 0.4em 0.5em 0.4em;
  box-sizing: border-box;
  padding-bottom: 0.4em;
  background-color: var(--sklearn-color-4);
}
#$id div.sk-label label {
  font-family: monospace;
  font-weight: bold;
  display: inline-block;
  line-height: 1.2em;
}
#$id div.sk-label-container {
  text-align: center;
}
#$id div.sk-container {
  /* jupyter's `normalize.less` sets `[hidden] { display: none; }`
     but bootstrap.min.css set `[hidden] { display: none !important; }`
     so we also need the `!important` here to be able to override the
     default hidden behavior on the sphinx rendered scikit-learn.org.
     See: https://github.com/scikit-learn/scikit-learn/issues/21755 */
  display: inline-block !important;
  position: relative;
}
#$id div.sk-text-repr-fallback {
  display: none;
}
""".replace("  ", "").replace("\n", "")  # noqa
