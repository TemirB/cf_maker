"""Render every page for inspection and check basic PDF integrity."""
from pathlib import Path
import json
import math
import pymupdf
from PIL import Image, ImageDraw

BASE = Path(__file__).resolve().parent
OUT = BASE / "build/review"
OUT.mkdir(parents=True, exist_ok=True)
doc = pymupdf.open(BASE / "presentation_19-10-2026_for_begginers.pdf")
metadata=json.loads((BASE/"build/talk_metadata.json").read_text())
assert len(doc) == metadata["total_slide_count"]
issues=[]
tiles=[]
for i,p in enumerate(doc):
    assert abs(p.rect.width/p.rect.height-16/9)<.001
    pix=p.get_pixmap(matrix=pymupdf.Matrix(2.5,2.5),alpha=False)
    pix.save(OUT/f"slide_{i+1:02d}.png")
    im=Image.frombytes("RGB",(pix.width,pix.height),pix.samples)
    im.thumbnail((600,338))
    tile=Image.new("RGB",(620,368),"#dce3e8")
    tile.paste(im,((620-im.width)//2,23))
    ImageDraw.Draw(tile).text((12,5),f"{i+1:02d}",fill="#142d43")
    tiles.append(tile)
    for block in p.get_text("dict")["blocks"]:
        if "lines" not in block: continue
        for line in block["lines"]:
            for span in line["spans"]:
                r=pymupdf.Rect(span["bbox"])
                if (r.x0 < -.2 or r.y0 < -.2 or r.x1>p.rect.width+.2 or r.y1>p.rect.height+.2):
                    issues.append({"page":i+1,"text":span["text"],"bbox":list(r)})
for start in range(0,len(tiles),9):
    sheet=Image.new("RGB",(3*620,3*368),"#dce3e8")
    for j,tile in enumerate(tiles[start:start+9]):sheet.paste(tile,((j%3)*620,(j//3)*368))
    sheet.save(OUT/f"contact_{start//9+1}.png")
contact=Image.new("RGB",(3*620,math.ceil(len(tiles)/3)*368),"#dce3e8")
for j,tile in enumerate(tiles):contact.paste(tile,((j%3)*620,(j//3)*368))
contact.save(BASE/"preview_contact.png")
(OUT/"integrity.json").write_text(json.dumps({"pages":len(doc),"aspect_ratio":"16:9","text_outside_page":issues},ensure_ascii=False,indent=2))
print(f"Rendered {len(doc)} slides; text outside page: {len(issues)}")
