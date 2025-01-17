import inspect
import os
import pydoc
import re
import textwrap

from jinja2 import FileSystemLoader
from jinja2.sandbox import SandboxedEnvironment
from sphinx.jinja2glue import BuiltinTemplateLoader

from .docscrape import ClassDoc, FunctionDoc, NumpyDocString, ObjDoc
from .docscrape import get_doc_object as get_doc_object_orig
from .xref import make_xref

IMPORT_MATPLOTLIB_RE = r"\b(import +matplotlib|from +matplotlib +import)\b"


class SphinxDocString(NumpyDocString):
    def __init__(self, docstring, config=None):
        if config is None:
            config = {}
        NumpyDocString.__init__(self, docstring, config=config)
        self.load_config(config)

    def load_config(self, config):
        self.use_plots = config.get("use_plots", False)
        self.class_members_toctree = config.get("class_members_toctree", True)
        self.attributes_as_param_list = config.get("attributes_as_param_list", True)
        self.xref_param_type = config.get("xref_param_type", False)
        self.xref_aliases = config.get("xref_aliases", dict())
        self.xref_ignore = config.get("xref_ignore", set())
        self.template = config.get("template", None)
        if self.template is None:
            template_dirs = [os.path.join(os.path.dirname(__file__), "templates")]
            template_loader = FileSystemLoader(template_dirs)
            template_env = SandboxedEnvironment(loader=template_loader)
            self.template = template_env.get_template("numpydoc_docstring.rst")

    # string conversion routines
    def _str_header(self, name, symbol="`"):
        return [".. rubric:: " + name, ""]

    def _str_field_list(self, name):
        return [":" + name + ":"]

    def _str_indent(self, doc, indent=4):
        out = []
        for line in doc:
            out += [" " * indent + line]
        return out

    def _str_signature(self):
        return [""]

    def _str_summary(self):
        return self["Summary"] + [""]

    def _str_extended_summary(self):
        return self["Extended Summary"] + [""]

    def _str_returns(self, name="Returns"):
        named_fmt = "**%s** : %s"
        unnamed_fmt = "%s"

        out = []
        if self[name]:
            out += self._str_field_list(name)
            out += [""]
            for param in self[name]:
                param_type = param.type
                if param_type and self.xref_param_type:
                    param_type = make_xref(
                        param_type, self.xref_aliases, self.xref_ignore
                    )
                if param.name:
                    out += self._str_indent(
                        [named_fmt % (param.name.strip(), param_type)]
                    )
                else:
                    out += self._str_indent([unnamed_fmt % param_type.strip()])
                if not param.desc:
                    out += self._str_indent([".."], 8)
                else:
                    out += self._str_indent(param.desc, 8)
                out += [""]
        return out

    def _escape_args_and_kwargs(self, name):
        if name[:2] == "**":
            return r"\*\*" + name[2:]
        elif name[:1] == "*":
            return r"\*" + name[1:]
        else:
            return name

    def _process_param(self, param, desc, fake_autosummary):
        """Determine how to display a parameter

        Emulates autosummary behavior if fake_autosummary

        Parameters
        ----------
        param : str
            The name of the parameter
        desc : list of str
            The parameter description as given in the docstring. This is
            ignored when autosummary logic applies.
        fake_autosummary : bool
            If True, autosummary-style behaviour will apply for params
            that are attributes of the class and have a docstring.

        Returns
        -------
        display_param : str
            The marked up parameter name for display. This may include a link
            to the corresponding attribute's own documentation.
        desc : list of str
            A list of description lines. This may be identical to the input
            ``desc``, if ``autosum is None`` or ``param`` is not a class
            attribute, or it will be a summary of the class attribute's
            docstring.

        Notes
        -----
        This does not have the autosummary functionality to display a method's
        signature, and hence is not used to format methods.  It may be
        complicated to incorporate autosummary's signature mangling, as it
        relies on Sphinx's plugin mechanism.
        """
        param = self._escape_args_and_kwargs(param.strip())
        # param = param.strip()
        # XXX: If changing the following, please check the rendering when param
        # ends with '_', e.g. 'word_'
        # See https://github.com/numpy/numpydoc/pull/144
        display_param = f"**{param}**"

        if not fake_autosummary:
            return display_param, desc

        param_obj = getattr(self._obj, param, None)
        if not (
            callable(param_obj)
            or isinstance(param_obj, property)
            or inspect.isgetsetdescriptor(param_obj)
            or inspect.ismemberdescriptor(param_obj)
        ):
            param_obj = None
        obj_doc = pydoc.getdoc(param_obj)

        if not (param_obj and obj_doc):
            return display_param, desc

        prefix = getattr(self, "_name", "")
        if prefix:
            link_prefix = f"{prefix}."
        else:
            link_prefix = ""

        # Referenced object has a docstring
        display_param = f":obj:`{param} <{link_prefix}{param}>`"
        if obj_doc:
            # Overwrite desc. Take summary logic of autosummary
            desc = re.split(r"\n\s*\n", obj_doc.strip(), 1)[0]
            # XXX: Should this have DOTALL?
            #      It does not in autosummary
            m = re.search(r"^([A-Z].*?\.)(?:\s|$)", " ".join(desc.split()))
            if m:
                desc = m.group(1).strip()
            else:
                desc = desc.partition("\n")[0]
            desc = desc.split("\n")
        return display_param, desc

    def _str_param_list(self, name, fake_autosummary=False):
        """Generate RST for a listing of parameters or similar

        Parameter names are displayed as bold text, and descriptions
        are in definition lists.

        Parameters
        ----------
        name : str
            Section name (e.g. Parameters)
        fake_autosummary : bool
            When True, the parameter names may correspond to attributes of the
            object beign documented, usually ``property`` instances on a class.
            In this case, names will be linked to fuller descriptions.

        Returns
        -------
        rst : list of str
        """
        out = []
        if self[name]:
            out += self._str_field_list(name)
            out += [""]
            for param in self[name]:
                display_param, desc = self._process_param(
                    param.name, param.desc, fake_autosummary
                )
                parts = []
                if display_param:
                    parts.append(display_param)
                param_type = param.type
                if param_type:
                    param_type = param.type
                    if self.xref_param_type:
                        param_type = make_xref(
                            param_type, self.xref_aliases, self.xref_ignore
                        )
                    parts.append(param_type)
                out += self._str_indent([" : ".join(parts)])

                if not desc:
                    # empty definition
                    desc = [".."]
                out += self._str_indent(desc, 8)
                out += [""]

        return out

    def _str_member_list(self, name):
        """
        Generate a member listing, autosummary:: table where possible,
        and a table where not.

        """
        out = []
        if self[name]:
            out += [f".. rubric:: {name}", ""]
            prefix = getattr(self, "_name", "")

            if prefix:
                prefix = f"~{prefix}."

            autosum = []
            others = []
            for param in self[name]:
                param = param._replace(name=param.name.strip())

                # Check if the referenced member can have a docstring or not
                param_obj = getattr(self._obj, param.name, None)
                if not (
                    callable(param_obj)
                    or isinstance(param_obj, property)
                    or inspect.isdatadescriptor(param_obj)
                ):
                    param_obj = None

                if param_obj and pydoc.getdoc(param_obj):
                    # Referenced object has a docstring
                    autosum += [f"   {prefix}{param.name}"]
                else:
                    others.append(param)

            if autosum:
                out += [".. autosummary::"]
                if self.class_members_toctree:
                    out += ["   :toctree:"]
                out += [""] + autosum

            if others:
                maxlen_0 = max(3, max(len(p.name) + 4 for p in others))
                hdr = "=" * maxlen_0 + "  " + "=" * 10
                fmt = "%%%ds  %%s  " % (maxlen_0,)
                out += ["", "", hdr]
                for param in others:
                    name = "**" + param.name.strip() + "**"
                    desc = " ".join(x.strip() for x in param.desc).strip()
                    if param.type:
                        desc = f"({param.type}) {desc}"
                    out += [fmt % (name, desc)]
                out += [hdr]
            out += [""]
        return out

    def _str_section(self, name):
        out = []
        if self[name]:
            out += self._str_header(name)
            content = textwrap.dedent("\n".join(self[name])).split("\n")
            out += content
            out += [""]
        return out

    def _str_see_also(self, func_role):
        out = []
        if self["See Also"]:
            see_also = super()._str_see_also(func_role)
            out = [".. seealso::", ""]
            out += self._str_indent(see_also[2:])
        return out

    def _str_warnings(self):
        out = []
        if self["Warnings"]:
            out = [".. warning::", ""]
            out += self._str_indent(self["Warnings"])
            out += [""]
        return out

    def _str_index(self):
        idx = self["index"]
        out = []
        if len(idx) == 0:
            return out

        out += [f".. index:: {idx.get('default', '')}"]
        for section, references in idx.items():
            if section == "default":
                continue
            elif section == "refguide":
                out += [f"   single: {', '.join(references)}"]
            else:
                out += [f"   {section}: {','.join(references)}"]
        out += [""]
        return out

    def _str_references(self):
        out = []
        if self["References"]:
            out += self._str_header("References")
            if isinstance(self["References"], str):
                self["References"] = [self["References"]]
            out.extend(self["References"])
            out += [""]
            # Latex collects all references to a separate bibliography,
            # so we need to insert links to it
            out += [".. only:: latex", ""]
            items = []
            for line in self["References"]:
                m = re.match(r".. \[([a-z0-9._-]+)\]", line, re.IGNORECASE)
                if m:
                    items.append(m.group(1))
            out += ["   " + ", ".join([f"[{item}]_" for item in items]), ""]
        return out

    def _str_examples(self):
        examples_str = "\n".join(self["Examples"])

        if (
            self.use_plots
            and re.search(IMPORT_MATPLOTLIB_RE, examples_str)
            and "plot::" not in examples_str
        ):
            out = []
            out += self._str_header("Examples")
            out += [".. plot::", ""]
            out += self._str_indent(self["Examples"])
            out += [""]
            return out
        else:
            return self._str_section("Examples")

    def __str__(self, indent=0, func_role="obj"):
        ns = {
            "signature": self._str_signature(),
            "index": self._str_index(),
            "summary": self._str_summary(),
            "extended_summary": self._str_extended_summary(),
            "parameters": self._str_param_list("Parameters"),
            "attributes": (
                self._str_param_list("Attributes", fake_autosummary=True)
                if self.attributes_as_param_list
                else self._str_member_list("Attributes")
            ),
            "methods": self._str_member_list("Methods"),
            "returns": self._str_returns("Returns"),
            "yields": self._str_returns("Yields"),
            "receives": self._str_returns("Receives"),
            "other_parameters": self._str_param_list("Other Parameters"),
            "raises": self._str_returns("Raises"),
            "warns": self._str_returns("Warns"),
            "warnings": self._str_warnings(),
            "see_also": self._str_see_also(func_role),
            "notes": self._str_section("Notes"),
            "references": self._str_references(),
            "examples": self._str_examples(),
        }
        ns = {k: "\n".join(v) for k, v in ns.items()}

        rendered = self.template.render(**ns)
        return "\n".join(self._str_indent(rendered.split("\n"), indent))


class SphinxFunctionDoc(SphinxDocString, FunctionDoc):
    def __init__(self, obj, doc=None, config=None):
        if config is None:
            config = {}
        self.load_config(config)
        FunctionDoc.__init__(self, obj, doc=doc, config=config)


class SphinxClassDoc(SphinxDocString, ClassDoc):
    def __init__(self, obj, doc=None, func_doc=None, config=None):
        if config is None:
            config = {}
        self.load_config(config)
        ClassDoc.__init__(self, obj, doc=doc, func_doc=None, config=config)


class SphinxObjDoc(SphinxDocString, ObjDoc):
    def __init__(self, obj, doc=None, config=None):
        if config is None:
            config = {}
        self.load_config(config)
        ObjDoc.__init__(self, obj, doc=doc, config=config)


def get_doc_object(obj, what=None, doc=None, config=None, builder=None):
    if config is None:
        config = {}

    template_dirs = [os.path.join(os.path.dirname(__file__), "templates")]
    if builder is not None:
        template_loader = BuiltinTemplateLoader()
        template_loader.init(builder, dirs=template_dirs)
    else:
        template_loader = FileSystemLoader(template_dirs)
    template_env = SandboxedEnvironment(loader=template_loader)
    config["template"] = template_env.get_template("numpydoc_docstring.rst")

    return get_doc_object_orig(
        obj,
        what=what,
        doc=doc,
        config=config,
        class_doc=SphinxClassDoc,
        func_doc=SphinxFunctionDoc,
        obj_doc=SphinxObjDoc,
    )
