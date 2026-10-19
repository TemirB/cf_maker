"""Build Beamer notes and an A4 speaker script from the editable Markdown source."""
import json
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
STEM = "19-10-2026_for_begginers"
raw = (BASE / f"speaker_notes_{STEM}.md").read_text()
parts = re.split(r"^## (\d{2})\. (.+)$", raw, flags=re.M)
notes = BASE / "notes"
notes.mkdir(exist_ok=True)

def escape(text):
    trans = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
             "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}"}
    return "".join(trans.get(c,c) for c in text)

MATH_TOKEN = re.compile(
    r"(?<![\\$])(\$\$.*?(?<!\\)\$\$|\$(?!\$)(?:\\.|[^$\\])+\$(?!\$))",
    flags=re.S,
)

def render_text(text):
    """Escape prose, retaining paragraphs while joining wrapped source lines."""
    text = re.sub(r"\n[ \t]*\n", "\n\n", text)
    text = re.sub(r"(?<!\n)\n(?!\n)", " ", text)
    chunks = re.split(r"(https?://[^\s)]+)", text)
    for index, chunk in enumerate(chunks):
        if chunk.startswith(("https://", "http://")):
            url = chunk.rstrip(".,;")
            chunks[index] = r"\url{" + url + "}" + escape(chunk[len(url):])
        else:
            chunks[index] = escape(chunk)
    text = "".join(chunks)
    return re.sub(r"(?m)^(Время|Тезис|Визуал|На экране|Текст|Пауза|Переход|Справка|Источники(?: и условия)?):",
                  r"\\textbf{\1:}", text)

def render_body(body):
    """Render $inline$ and $$display$$ math; keep all other input as prose."""
    rendered = []
    position = 0
    for match in MATH_TOKEN.finditer(body):
        rendered.append(render_text(body[position:match.start()]))
        formula = match.group(0)
        if formula.startswith("$$"):
            rendered.append(r"\[" + formula[2:-2] + r"\]")
        else:
            rendered.append(r"\(" + formula[1:-1] + r"\)")
        position = match.end()
    rendered.append(render_text(body[position:]))
    return "".join(rendered)

sections=[]
seconds=0
main_slides=[]
reserve_slides=[]
for i in range(1,len(parts),3):
    number,title,body=parts[i:i+3]
    content=render_body(body.strip())+"\n"
    (notes/f"n{number}.tex").write_text(content)
    page_break = r"\clearpage"+"\n" if sections else ""
    sections.append(page_break+r"\section*{Слайд "+number+". "+escape(title)+"}\n"+content)
    tm=re.search(r"Время: (\d+):(\d+)",body)
    if tm:
        assert not reserve_slides, "Timed slides must precede the reserve"
        duration = int(tm[1])*60+int(tm[2])
        seconds += duration
        main_slides.append({"number": int(number), "title": title, "seconds": duration})
    else:
        reserve_slides.append(int(number))

header=r"""\documentclass[11pt,a4paper]{article}
\usepackage{fontspec,polyglossia}
\usepackage{amsmath,amssymb}
\newcommand{\hc}{\hbar c}
\newcommand{\dd}{\mathrm d}
\setdefaultlanguage{russian}
\setmainfont{Arial}
\newfontfamily\cyrillicfont{Arial}[Script=Cyrillic]
\newfontfamily\cyrillicfonttt{Arial}[Script=Cyrillic]
\usepackage[margin=21mm]{geometry}
\usepackage{xcolor,hyperref}
\definecolor{ink}{HTML}{142D43}
\definecolor{teal}{HTML}{007F86}
\hypersetup{colorlinks=true,urlcolor=teal,pdftitle={Заметки докладчика: фемтоскопия},pdfauthor={Балинов Темир}}
\urlstyle{same}
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
\linespread{1.08}
\setlength{\emergencystretch}{3em}
\begin{document}
{\LARGE\bfseries\color{ink} От импульсов к пространству}\par
{\large Заметки докладчика}\par
Балинов Темир, физический факультет МГУ. 19.10.2026.\par
__TALK_SUMMARY__\par
\medskip
"""
assert main_slides and reserve_slides, "Expected a main talk and reserve"
assert [int(parts[i]) for i in range(1,len(parts),3)] == list(range(1,len(sections)+1))
assert 40*60 <= seconds <= 45*60, (seconds, "Talk must fit 40–45 minutes")
summary = (f"{len(main_slides)} основных слайдов: {seconds//60}:{seconds%60:02d} "
           f"с объяснением формул и паузами. Резерв: слайды {reserve_slides[0]}–{reserve_slides[-1]}.")
header = header.replace("__TALK_SUMMARY__", summary)
(BASE/f"speaker_notes_{STEM}.tex").write_text(header+"\n".join(sections)+"\n\\end{document}\n")
(BASE/"build").mkdir(exist_ok=True)
(BASE/"build/talk_metadata.tex").write_text(f"\\newcommand{{\\MainSlideCount}}{{{len(main_slides)}}}\n")
(BASE/"build/talk_metadata.json").write_text(json.dumps({
    "main_slide_count":len(main_slides), "total_slide_count":len(sections),
    "duration_seconds":seconds, "slides":main_slides, "reserve":reserve_slides,
}, ensure_ascii=False, indent=2)+"\n")
print(f"Prepared {len(sections)} slide notes; main talk: {seconds//60}:{seconds%60:02d}.")
