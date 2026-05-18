#모듈 불러오기
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter, NullFormatter # <-- 이 줄이 빠져서 발생하는 에러입니다.
import streamlit as st
import pandas as pd
import io
import os
import ssl
import subprocess
import requests
import matplotlib.pyplot as plt
import numpy as np
import json
from Bio.Seq import Seq  # siRNA 역상보 서열 계산용
from Bio.SeqUtils import MeltingTemp as mt
import free_energy
import general_helpers

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML

from stmol import showmol
import py3Dmol


#보안 및 기본 설정
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "csw1567@naver.com" 

if 'setup_done' not in st.session_state:
    st.set_page_config(
        page_title="RNAi design prototype(KNU)",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    st.session_state['setup_done'] = True




#모든 텝에 적용될 CSS 디자인 넣기
st.markdown("""
    <style>
    /* 1. 탭 전체 글자 크기 및 스타일 */
    .stTabs [data-baseweb="tab-list"] button p {
        font-size: 18px;
        font-weight: bold;
    }

    /* 2. 탭 버튼 사이에 세로선 넣기 */
    .stTabs [data-baseweb="tab-list"] button {
        border-right: 1px solid #d3d3d3; /* 오른쪽에 회색 세로선 추가 */
        border-radius: 0px; /* 각진 모양으로 변경하여 선이 깔끔하게 보이게 함 */
        padding-left: 20px;
        padding-right: 20px;
    }

    /* 3. 마지막 탭은 세로선 없애기 */
    .stTabs [data-baseweb="tab-list"] button:last-child {
        border-right: none;
    }

    /* 4. 선택된 탭 배경색 및 강조 */
    .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
        background-color: #e6f3ff;
        color: #1a2a6c;
    }
    </style>
    """, unsafe_allow_html=True)

# 4. 사이드바 구성 
with st.sidebar:
    st.title("KNU Bio-Lab")
    st.caption("Insect Physiology Laboratory")
    st.markdown("---")

    # 문의처 섹션 (HTML/CSS 사용)
    st.markdown("""
        <div style="background-color:#f0f2f6; padding:15px; border-radius:10px; border-left:5px solid #1a2a6c; margin-bottom:20px;">
            <div style="color:#1a2a6c; font-weight:bold; margin-bottom:5px;">📧 Inquiry & Support</div>
            <p style="margin:0; font-size:0.95rem;"><b>곤충생리학연구실</b></p>
            <p style="margin:0; font-size:0.9rem;">최성환 연구원에게 찾아오세요.</p>
        </div>
    """, unsafe_allow_html=True)
    
    # 유튜브 영상 섹션 
    st.markdown("### 📺 RNAi Learning Resource")
    st.video("https://www.youtube.com/watch?v=zVe3Nom2SpI")
    st.caption("Mechanism of RNA Interference")
    st.markdown("---")

    
    
    # 버전 정보 (맨 밑)
    st.markdown("""
        <div style="font-size:0.8rem; color:grey; text-align:center;">
            System Version: 1.1.0-beta<br>
            © 2026 KNU Insect Physiology Lab
        </div>
    """, unsafe_allow_html=True)


#메인 타이틀 디자인
st.title("Bursaphelenchus xylophilus genetic analysis tool")
tab1, tab2, tab3, tab4, tab5, = st.tabs([
    "타겟 유전자 정보 조회",
    "RNAi디자인 및 off-target",
    "가상 gel imagine 생성",
    "농도 계산기",
    "Gene조정 및 파일형식 변환"
])


import re
import urllib.parse
import subprocess
import tempfile
import re






with tab1:
    st.header("프라이머를 통한 타겟 유전자 분석")

    # ──────────────────────────────────────────────
    # 헬퍼 함수
    # ──────────────────────────────────────────────

    @st.cache_data
    def build_id_mapping_table():
        current_dir  = os.path.dirname(os.path.abspath(__file__))
        protein_path = os.path.join(current_dir, "pwn_pro_named.fa")
        mapping = {}
        if os.path.exists(protein_path):
            for record in SeqIO.parse(protein_path, "fasta"):
                prot_id    = str(record.id)
                desc_parts = record.description.split(" ")
                if len(desc_parts) >= 3 and desc_parts[1] == prot_id:
                    name = " ".join(desc_parts[2:]).strip()
                elif len(desc_parts) >= 2:
                    name = " ".join(desc_parts[1:]).strip()
                else:
                    name = "Hypothetical Protein"
                name = name.replace("%0A", " ").replace("%0a", " ").strip()
                if len(name) > 60:
                    name = name[:60] + "..."
                mapping[prot_id]               = name
                mapping[prot_id.split(".")[0]] = name
        return mapping

    @st.cache_resource
    def ensure_blast_db():
        current_dir = os.path.dirname(os.path.abspath(__file__))
        db_path     = os.path.join(current_dir, "pwn_blast_db", "pwn_blast_db")
        return db_path

    def get_protein_name(locus_id, mapping):
        raw  = str(locus_id).strip()
        base = raw.split(".")[0]
        return mapping.get(raw) or mapping.get(base) or "Hypothetical Protein"

    WBPS_SPECIES = "bursaphelenchus_xylophilus_prjea64437"

    # ──────────────────────────────────────────────
    # 섹션 1: 프라이머 BLAST 분석
    # ──────────────────────────────────────────────
    st.subheader("1. 프라이머 서열 입력")

    query_seq = st.text_area(
        "프라이머 서열을 입력하세요 (ATGC...)",
        height=100,
        placeholder="예: GCTTCACAGCAGTGCCAAG",
        key="tab1_primer_input"
    )

    run_btn = st.button("타겟 유전자 분석 실행", use_container_width=True, type="primary")

    if run_btn:
        if not query_seq.strip():
            st.warning("프라이머 서열을 입력해 주세요.")
        else:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            temp_query  = os.path.join(current_dir, "temp_query.fa")
            result_csv  = os.path.join(current_dir, "blast_result.csv")

            with st.spinner("BLAST 분석 중..."):
                try:
                    db_path = ensure_blast_db()

                    with open(temp_query, "w") as f:
                        f.write(f">Query\n{query_seq.strip()}")

                    subprocess.run([
                        "blastn", "-query", temp_query, "-db", db_path,
                        "-out", result_csv,
                        "-outfmt", "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                        "-task", "blastn-short", "-evalue", "1000"
                    ], check=True)

                    if os.path.exists(result_csv) and os.path.getsize(result_csv) > 0:
                        df = pd.read_csv(result_csv, names=[
                            "Query", "Locus ID", "Identity(%)", "Length",
                            "Mismatch", "Gaps", "Q_Start", "Q_End",
                            "S_Start", "S_End", "E-value", "BitScore"
                        ])
                        mapping = build_id_mapping_table()
                        df["Protein Name"] = df["Locus ID"].apply(
                            lambda x: get_protein_name(x, mapping)
                        )
                        df_sorted = df.sort_values("E-value").reset_index(drop=True)
                        df_sorted.index += 1
                        st.session_state["blast_df"]   = df_sorted
                        st.session_state["blast_done"] = True
                    else:
                        st.error("BLAST 결과가 없습니다. 서열을 확인해 주세요.")
                        st.session_state["blast_done"] = False

                except Exception as e:
                    st.error(f"분석 중 오류 발생: {e}")
                    st.session_state["blast_done"] = False

    # ──────────────────────────────────────────────
    # 섹션 2: BLAST 결과 테이블
    # ──────────────────────────────────────────────
    if st.session_state.get("blast_done"):
        df_sorted = st.session_state["blast_df"]

        st.markdown("---")
        st.subheader("2. 분석 결과 — 후보 유전자 목록")
        st.caption(f"총 {len(df_sorted)}개 결과 | E-value 오름차순 정렬 (상위 = 높은 유사도)")

        st.dataframe(
            df_sorted[["Protein Name", "Locus ID", "Identity(%)", "E-value", "BitScore"]],
            use_container_width=True,
            height=min(400, 40 + len(df_sorted) * 35)
        )

        csv_data = df_sorted[["Protein Name", "Locus ID", "Identity(%)", "E-value", "BitScore"]].to_csv(index=True)
        st.download_button(
            "결과 CSV 다운로드",
            data=csv_data,
            file_name="blast_result_named.csv",
            mime="text/csv"
        )

    # ──────────────────────────────────────────────
    # 섹션 3: WormBase ParaSite 조회
    # ──────────────────────────────────────────────
    st.markdown("---")
    st.subheader("3. WormBase ParaSite 조회")
    st.info("BLAST 결과에서 선택하거나 Locus ID를 직접 입력하세요.")

    col_select, col_manual = st.columns([2, 1])

    with col_select:
        if st.session_state.get("blast_done"):
            df_ref  = st.session_state["blast_df"]
            options = ["— 직접 입력 —"] + [
                f"{row['Protein Name']}  [{row['Locus ID']}]"
                for _, row in df_ref.iterrows()
            ]
            selected = st.selectbox("결과에서 선택", options)
            if selected != "— 직접 입력 —":
                match   = re.search(r'\[(.+?)\]', selected)
                auto_id = match.group(1) if match else ""
            else:
                auto_id = ""
        else:
            st.caption("위에서 BLAST 분석을 먼저 실행하거나 아래에 직접 입력하세요.")
            auto_id = ""

    with col_manual:
        manual_id = st.text_input(
            "직접 입력 (우선 적용)",
            placeholder="예: BXY_0416800.1",
            key="wbps_manual_id"
        )

    query_id = manual_id.strip() if manual_id.strip() else auto_id
    base_id  = query_id.split(".")[0] if query_id else ""

    if query_id:
        st.caption(f"조회 대상: `{base_id}`")

    link_btn = st.button("WormBase ParaSite에서 검색", use_container_width=True)

    if link_btn:
        if not base_id:
            st.warning("조회할 ID를 선택하거나 입력해 주세요.")
        else:
            mapping   = build_id_mapping_table()
            prot_name = get_protein_name(query_id, mapping)
            st.markdown(
                f"[**{base_id}** — WormBase ParaSite 전체 검색]"
                f"(https://parasite.wormbase.org/Multi/Search/Results?q={base_id})"
            )
            st.caption(f"단백질 이름: {prot_name}")
                
with tab2:
    st.header("siRNA 디자인 및 Off-target 분석")

    # ──────────────────────────────────────────────
    # 공통 함수
    # ──────────────────────────────────────────────

    def generate_sirnas(sequence, sirna_size=21):
        """서열에서 siRNA 후보군 생성"""
        sirnas   = []
        sequence = sequence.upper().replace("U", "T")
        for i in range(len(sequence) - sirna_size + 1):
            sirna_seq = sequence[i:i + sirna_size]
            if len(sirna_seq) == sirna_size:
                sirnas.append({
                    "position": i + 1,
                    "sequence": sirna_seq,
                    "name":     f"siRNA_{i+1}"
                })
        return sirnas

    def calc_gc_content(seq):
        """GC 함량 계산"""
        seq = seq.upper()
        gc  = seq.count('G') + seq.count('C')
        return round(gc / len(seq) * 100, 1) if seq else 0

    def calc_free_energy(seq):
        """말단 자유에너지 근사 계산 (nearest-neighbor)"""
        nn_params = {
            'AA': -1.0, 'AT': -0.9, 'TA': -0.6, 'CA': -1.7,
            'GT': -1.5, 'CT': -1.3, 'GA': -1.6, 'CG': -3.6,
            'GC': -3.1, 'GG': -3.1, 'AC': -1.8, 'TC': -1.3,
            'AG': -1.5, 'TG': -1.7, 'TT': -1.0, 'CC': -3.1
        }
        energy = 0.0
        for i in range(len(seq) - 1):
            pair    = seq[i:i+2].upper()
            energy += nn_params.get(pair, -1.5)
        return round(energy, 2)

    def check_strand_selection(sense_seq, end_nucleotides=3):
        """가닥 선택 규칙 체크"""
        sense_5prime     = sense_seq[:end_nucleotides]
        antisense_5prime = sense_seq[-end_nucleotides:]
        sense_energy     = calc_free_energy(sense_5prime)
        antisense_energy = calc_free_energy(antisense_5prime)
        return antisense_energy <= sense_energy, sense_energy, antisense_energy

    def check_gc_rule(seq):
        """GC 함량 30~70% 체크"""
        gc = calc_gc_content(seq)
        return 30 <= gc <= 70, gc

    def check_terminal_rule(seq):
        """말단 뉴클레오타이드 규칙 체크"""
        pos1  = seq[0]  in ['A', 'T']
        pos19 = seq[-3] in ['A', 'T']
        return pos1 and pos19

    def check_homopolymer(seq, max_run=4):
        """호모폴리머 체크"""
        for base in ['A', 'T', 'G', 'C']:
            if base * max_run in seq:
                return False
        return True

    def calculate_efficiency(sirna_seq):
        """종합 효율 점수 계산 (0~100)"""
        score  = 0
        detail = {}

        gc_ok, gc = check_gc_rule(sirna_seq)
        if gc_ok:
            score += 30
        detail["GC 함량"] = f"{gc}% {'통과' if gc_ok else '실패'}"

        strand_ok, s_e, as_e = check_strand_selection(sirna_seq)
        if strand_ok:
            score += 30
        detail["가닥 선택"] = f"Sense {s_e} / Antisense {as_e} kcal/mol {'통과' if strand_ok else '실패'}"

        terminal_ok = check_terminal_rule(sirna_seq)
        if terminal_ok:
            score += 20
        detail["말단 규칙"] = f"{'통과' if terminal_ok else '실패'}"

        homo_ok = check_homopolymer(sirna_seq)
        if homo_ok:
            score += 20
        detail["호모폴리머"] = f"{'없음' if homo_ok else '있음'}"

        return score, detail

    def get_antisense(seq):
        """siRNA antisense 서열 (역상보)"""
        return str(Seq(seq).reverse_complement())

    def run_blast(query_seq, db_path, work_dir, evalue="10", word_size="7"):
        """BLAST로 off-target 검색"""
        tmp_fa  = os.path.join(work_dir, "sirna_query_tmp.fa")
        tmp_out = os.path.join(work_dir, "sirna_blast_tmp.csv")
        with open(tmp_fa, "w") as f:
            f.write(f">query\n{query_seq}\n")
        subprocess.run([
            "blastn",
            "-query",     tmp_fa,
            "-db",        db_path,
            "-out",       tmp_out,
            "-outfmt",    "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-task",      "blastn-short",
            "-evalue",    evalue,
            "-word_size", word_size
        ], check=True)
        if os.path.exists(tmp_out) and os.path.getsize(tmp_out) > 0:
            df = pd.read_csv(tmp_out, names=[
                "Query", "Target", "Identity(%)", "Length",
                "Mismatch", "Gap", "Q_Start", "Q_End",
                "S_Start", "S_End", "E-value", "BitScore"
            ])
            return df
        return pd.DataFrame()

    def parse_fasta_or_plain(text):
        """FASTA 또는 순수 서열 파싱"""
        lines   = text.strip().split("\n")
        results = []
        if any(l.startswith(">") for l in lines):
            name, seq = "", ""
            for l in lines:
                if l.startswith(">"):
                    if seq:
                        results.append({"name": name, "sequence": seq.upper()})
                    name, seq = l[1:].strip(), ""
                else:
                    seq += l.strip()
            if seq:
                results.append({"name": name, "sequence": seq.upper()})
        else:
            for i, l in enumerate(lines):
                l = l.strip()
                if l:
                    results.append({"name": f"siRNA_{i+1}", "sequence": l.upper()})
        return results

    # ──────────────────────────────────────────────
    # DB 자동 감지 (공통)
    # ──────────────────────────────────────────────
    current_dir = os.path.dirname(os.path.abspath(__file__))

    db_candidates = {
        "소나무재선충 (주타겟)": "pwn_blast_db/pwn_blast_db",
        "꿀벌 (Apis mellifera)": "honeybee_db/honeybee_db",
        "소나무 (Pinus)":        "pine_db/pine_db",
        "사람 (Homo sapiens)":   "human_db/human_db",
        "C. elegans":            "celegans_db/celegans_db",
    }

    available_dbs   = {}
    unavailable_dbs = []
    for name, path in db_candidates.items():
        full_path = os.path.join(current_dir, path)
        if os.path.exists(full_path + ".nhr"):
            available_dbs[name] = full_path
        else:
            unavailable_dbs.append(name)

    # ──────────────────────────────────────────────
    # 모드 선택
    # ──────────────────────────────────────────────
    mode = st.radio(
        "분석 모드 선택",
        ["RNAi Design Mode", "Off-target Prediction Mode"],
        horizontal=True
    )
    st.markdown("---")

    # ══════════════════════════════════════════════
    # MODE 0: RNAi Design Mode
    # ══════════════════════════════════════════════
    if "RNAi Design" in mode:
        st.subheader("RNAi Design Mode")
        st.info("mRNA 또는 CDS 서열을 입력하면 최적 siRNA 후보군을 설계합니다.")

        # 1. 서열 입력
        st.markdown("### 1. 서열 입력")
        uploaded = st.file_uploader(
            "CDS/mRNA FASTA 파일 업로드",
            type=["fasta", "fa", "fna", "txt"],
            key="design_upload"
        )
        target_fasta = ""
        if uploaded:
            target_fasta = uploaded.read().decode("utf-8")
            st.success(f"업로드 완료: {uploaded.name}")
        else:
            target_fasta = st.text_area(
                "또는 서열 직접 입력 (FASTA 또는 순수 서열)",
                height=150,
                placeholder=">gene_name\nATGCGTATCGATCG...",
                key="design_seq_input"
            )

        # 2. 파라미터
        st.markdown("### 2. 파라미터 설정")
        col_p1, col_p2 = st.columns(2)
        with col_p1:
            sirna_size = st.number_input(
                "siRNA 크기 (nt)",
                value=21, min_value=19, max_value=25,
                help="일반적으로 21nt 사용. 19~25nt 범위에서 조절 가능."
            )
        with col_p2:
            min_score = st.slider(
                "최소 효율 점수 필터",
                0, 100, 50, step=10,
                help="0~100점 중 이 점수 이상인 siRNA만 표시됩니다.\n"
                     "100점: 모든 규칙 통과\n"
                     "50점: 절반 이상 통과 (권장)\n"
                     "0점: 전체 후보 표시"
            )

        st.markdown("---")

        # 3. 실행
        design_btn = st.button("siRNA 설계 시작", use_container_width=True, type="primary")

        if design_btn:
            if not target_fasta.strip():
                st.error("서열을 입력하거나 파일을 업로드해 주세요.")
            else:
                raw = target_fasta.strip()
                if raw.startswith(">"):
                    lines    = raw.split("\n")
                    seq_name = lines[0][1:].strip()
                    sequence = "".join(lines[1:]).upper().replace(" ", "")
                else:
                    seq_name = "Query"
                    sequence = raw.upper().replace(" ", "")
                sequence = "".join(sequence.split())

                if len(sequence) < sirna_size:
                    st.error(f"서열이 너무 짧습니다. 최소 {sirna_size}nt 이상 입력해 주세요.")
                else:
                    with st.spinner("siRNA 후보군 분석 중..."):
                        sirnas  = generate_sirnas(sequence, sirna_size)
                        results = []
                        for s in sirnas:
                            score, detail = calculate_efficiency(s["sequence"])
                            gc = calc_gc_content(s["sequence"])
                            _, s_e, as_e = check_strand_selection(s["sequence"])
                            results.append({
                                "위치":         s["position"],
                                "siRNA 서열":   s["sequence"],
                                "Antisense":    get_antisense(s["sequence"]),
                                "효율 점수":    score,
                                "GC 함량(%)":   gc,
                                "Sense dG":     s_e,
                                "Antisense dG": as_e,
                                "GC 규칙":      "통과" if 30 <= gc <= 70 else "실패",
                                "가닥 선택":    "통과" if check_strand_selection(s["sequence"])[0] else "실패",
                                "말단 규칙":    "통과" if check_terminal_rule(s["sequence"]) else "실패",
                                "호모폴리머":   "없음" if check_homopolymer(s["sequence"]) else "있음",
                                "_detail":      detail
                            })

                        df = pd.DataFrame(results)
                        df = df[df["효율 점수"] >= min_score].sort_values(
                            "효율 점수", ascending=False
                        ).reset_index(drop=True)
                        df.index += 1
                        st.session_state["design_df"]   = df
                        st.session_state["design_done"] = True

        # 4. 결과 출력
        if st.session_state.get("design_done"):
            df = st.session_state["design_df"]

            st.markdown("---")
            st.subheader("siRNA 후보군 결과")
            st.caption(
                f"효율 점수 {min_score}점 이상 후보: {len(df)}개 | "
                f"효율 점수 내림차순 정렬"
            )

            score_counts = df["효율 점수"].value_counts().sort_index()
            st.bar_chart(score_counts)

            display_cols = ["위치", "siRNA 서열", "효율 점수", "GC 함량(%)",
                            "GC 규칙", "가닥 선택", "말단 규칙", "호모폴리머"]
            st.dataframe(df[display_cols], use_container_width=True)

            st.download_button(
                "결과 CSV 다운로드",
                data=df[["위치", "siRNA 서열", "Antisense", "효율 점수",
                          "GC 함량(%)", "Sense dG", "Antisense dG"]].to_csv(index=True),
                file_name="sirna_design.csv",
                mime="text/csv"
            )

            # 5. 개별 상세 및 off-target
            st.markdown("---")
            st.subheader("후보 siRNA 상세 및 Off-target 확인")

            if not df.empty:
                options = [
                    f"위치 {r['위치']} --- {r['siRNA 서열']} (점수: {r['효율 점수']})"
                    for _, r in df.iterrows()
                ]
                sel = st.selectbox("siRNA 선택", options, key="design_sel")
                idx = options.index(sel)
                row = df.iloc[idx]

                col_d1, col_d2 = st.columns(2)
                with col_d1:
                    st.markdown("**Sense (5->3):**")
                    st.code(row["siRNA 서열"], language="text")
                    st.markdown("**Antisense (3->5):**")
                    st.code(row["Antisense"], language="text")
                with col_d2:
                    st.metric("효율 점수", f"{row['효율 점수']} / 100")
                    for k, v in row["_detail"].items():
                        st.markdown(f"- **{k}:** {v}")

                st.markdown("#### Off-target 검색 파라미터")
                col_o1, col_o2, col_o3 = st.columns(3)
                with col_o1:
                    design_mismatch = st.slider(
                        "위험 기준 미스매치 수",
                        0, 5, 2, key="design_mismatch",
                        help="이 수 이하 미스매치를 위험으로 분류"
                    )
                with col_o2:
                    design_evalue = st.selectbox(
                        "E-value", ["1", "10", "100"], index=1,
                        key="design_evalue"
                    )
                with col_o3:
                    design_wordsize = st.selectbox(
                        "Word Size", ["7", "8", "11"], index=0,
                        key="design_wordsize"
                    )

                if available_dbs:
                    design_selected_dbs = st.multiselect(
                        "검색할 DB 선택",
                        list(available_dbs.keys()),
                        default=list(available_dbs.keys())[:1],
                        key="design_db_select"
                    )
                else:
                    st.warning("사용 가능한 DB가 없습니다.")
                    design_selected_dbs = []

                offtarget_btn = st.button(
                    "Off-target 확인",
                    use_container_width=True,
                    key="design_offtarget_btn"
                )

                if offtarget_btn and design_selected_dbs:
                    for db_name in design_selected_dbs:
                        db_path = available_dbs[db_name]
                        st.markdown(f"#### {db_name} 결과")
                        with st.spinner(f"{db_name} 검색 중..."):
                            try:
                                df_off = run_blast(
                                    row["siRNA 서열"], db_path, current_dir,
                                    evalue=design_evalue, word_size=design_wordsize
                                )
                                if not df_off.empty:
                                    if "재선충" in db_name:
                                        mapping = build_id_mapping_table()
                                        df_off["Protein Name"] = df_off["Target"].apply(
                                            lambda x: get_protein_name(x, mapping)
                                        )
                                    else:
                                        df_off["Protein Name"] = "-"

                                    df_off["위험도"] = df_off["Mismatch"].apply(
                                        lambda x: "위험" if int(x) <= design_mismatch else "낮음"
                                    )
                                    danger = len(df_off[df_off["위험도"] == "위험"])

                                    if danger > 0:
                                        st.warning(f"위험 off-target {danger}개 / 전체 {len(df_off)}개")
                                    else:
                                        st.success(f"위험 off-target 없음 / 전체 {len(df_off)}개")

                                    st.dataframe(
                                        df_off[["Protein Name", "Target", "Identity(%)",
                                                "Mismatch", "E-value", "위험도"]].sort_values("Mismatch"),
                                        use_container_width=True
                                    )
                                else:
                                    st.success("Off-target 없음")
                            except Exception as e:
                                st.error(f"검색 실패: {e}")

    # ══════════════════════════════════════════════
    # MODE 1: Off-target Prediction Mode
    # ══════════════════════════════════════════════
    else:
        st.subheader("Off-target Prediction Mode")
        st.info("보유한 siRNA 서열을 입력하고 검색할 DB를 선택하면 off-target을 예측합니다.")

        # 1. siRNA 서열 입력
        st.markdown("### 1. siRNA 서열 입력")
        sirna_input = st.text_area(
            "siRNA 서열 입력 (여러 개는 줄바꿈으로 구분, FASTA 형식도 가능)",
            height=150,
            placeholder="순수 서열:\nGGGATGGCTCAAAGGCGTAGT\nATCGATCGATCGATCGATCGA\n\nFASTA 형식:\n>siRNA_1\nGGGATGGCTCAAAGGCGTAGT",
            key="offtarget_input"
        )

        # 2. BLAST 파라미터
        st.markdown("### 2. BLAST 파라미터 설정")
        col_b1, col_b2, col_b3 = st.columns(3)
        with col_b1:
            mismatches_allowed = st.slider(
                "위험 off-target 기준 미스매치 수",
                min_value=0, max_value=5, value=2,
                help="이 수 이하의 미스매치를 위험 off-target으로 분류합니다.\n"
                     "0: 완벽 일치만 위험\n"
                     "2: 2개 이하 미스매치 위험 (권장)\n"
                     "5: 5개 이하 미스매치 위험 (엄격)"
            )
        with col_b2:
            evalue_threshold = st.selectbox(
                "E-value 임계값",
                options=["1", "10", "100", "1000"],
                index=1,
                help="높을수록 더 많은 히트 포함\n10: 권장 (off-target용)"
            )
        with col_b3:
            word_size = st.selectbox(
                "Word Size (검색 민감도)",
                options=["7", "8", "11"],
                index=0,
                help="작을수록 더 민감하게 검색\n7: 권장 (21nt siRNA용)"
            )

        # 3. DB 선택
        st.markdown("### 3. 검색할 DB 선택")
        if available_dbs:
            selected_dbs = st.multiselect(
                "검색할 DB 선택 (다중 선택 가능)",
                list(available_dbs.keys()),
                default=list(available_dbs.keys()),
                key="offtarget_db_select"
            )
        else:
            st.error("사용 가능한 DB가 없습니다.")
            selected_dbs = []

        if unavailable_dbs:
            with st.expander("DB 파일 없음 --- 추가 방법 보기"):
                for name in unavailable_dbs:
                    st.markdown(f"- {name}")
                st.code(
                    "makeblastdb -in honeybee_cds.fna -dbtype nucl -out honeybee_db/honeybee_db\n"
                    "makeblastdb -in human_cds.fna    -dbtype nucl -out human_db/human_db\n"
                    "makeblastdb -in pine_cds.fna     -dbtype nucl -out pine_db/pine_db",
                    language="bash"
                )

        # 4. 커스텀 DB
        st.markdown("### 4. 커스텀 DB 추가 (선택사항)")
        col_c1, col_c2 = st.columns([2, 1])
        with col_c1:
            custom_db_path = st.text_input(
                "커스텀 DB 경로",
                placeholder="예: custom_db/myorganism_db"
            )
        with col_c2:
            custom_db_name = st.text_input(
                "DB 이름",
                placeholder="예: 나방 (Bombyx mori)"
            )
        if custom_db_path and custom_db_name:
            full_custom = os.path.join(current_dir, custom_db_path)
            if os.path.exists(full_custom + ".nhr"):
                available_dbs[f"{custom_db_name}"] = full_custom
                st.success(f"커스텀 DB 추가됨: {custom_db_name}")
            else:
                st.error("해당 경로에 DB 파일이 없습니다.")

        st.markdown("---")

        # 5. 실행
        offtarget_run_btn = st.button(
            "Off-target 예측 시작",
            use_container_width=True,
            type="primary"
        )

        if offtarget_run_btn:
            if not sirna_input.strip():
                st.error("siRNA 서열을 입력해 주세요.")
            elif not selected_dbs:
                st.error("검색할 DB를 선택해 주세요.")
            else:
                sirnas_to_check = parse_fasta_or_plain(sirna_input)
                summary_all     = []

                st.markdown("---")
                st.subheader("분석 결과")

                for s in sirnas_to_check:
                    st.markdown(f"## {s['name']} --- {s['sequence']}")

                    score, detail = calculate_efficiency(s["sequence"])
                    gc = calc_gc_content(s["sequence"])
                    _, s_e, as_e = check_strand_selection(s["sequence"])

                    with st.expander("siRNA 효율 정보", expanded=True):
                        col_e1, col_e2, col_e3, col_e4 = st.columns(4)
                        with col_e1:
                            st.metric("효율 점수", f"{score} / 100")
                        with col_e2:
                            st.metric("GC 함량", f"{gc}%")
                        with col_e3:
                            st.metric("Sense dG", f"{s_e} kcal/mol")
                        with col_e4:
                            st.metric("Antisense dG", f"{as_e} kcal/mol")

                        col_r1, col_r2, col_r3, col_r4 = st.columns(4)
                        with col_r1:
                            st.markdown(f"GC 규칙: {detail['GC 함량']}")
                        with col_r2:
                            st.markdown(f"가닥 선택: {detail['가닥 선택']}")
                        with col_r3:
                            st.markdown(f"말단 규칙: {detail['말단 규칙']}")
                        with col_r4:
                            st.markdown(f"호모폴리머: {detail['호모폴리머']}")

                    st.markdown(f"**Antisense (3->5):** `{get_antisense(s['sequence'])}`")

                    sirna_summary = {
                        "siRNA":      s["name"],
                        "서열":       s["sequence"],
                        "효율 점수":  score,
                        "GC 함량(%)": gc,
                    }
                    total_danger    = 0
                    total_offtarget = 0

                    for db_name in selected_dbs:
                        db_path = available_dbs[db_name]
                        st.markdown(f"#### {db_name}")
                        try:
                            df_off = run_blast(
                                s["sequence"], db_path, current_dir,
                                evalue=evalue_threshold, word_size=word_size
                            )
                            if not df_off.empty:
                                if "재선충" in db_name:
                                    mapping = build_id_mapping_table()
                                    df_off["Protein Name"] = df_off["Target"].apply(
                                        lambda x: get_protein_name(x, mapping)
                                    )
                                else:
                                    df_off["Protein Name"] = "-"

                                df_off["위험도"] = df_off["Mismatch"].apply(
                                    lambda x: "위험" if int(x) <= mismatches_allowed else "낮음"
                                )
                                df_off["Seed 위험"] = df_off.apply(
                                    lambda row: "Seed 위험" if (
                                        int(row["Mismatch"]) <= mismatches_allowed and
                                        int(row["Q_Start"]) <= 2 and
                                        int(row["Q_End"]) >= 8
                                    ) else "-", axis=1
                                )

                                danger_df    = df_off[df_off["위험도"] == "위험"]
                                danger_count = len(danger_df)
                                total_danger    += danger_count
                                total_offtarget += len(df_off)

                                if danger_count > 0:
                                    st.warning(f"위험 off-target {danger_count}개 / 전체 {len(df_off)}개")
                                else:
                                    st.success(f"위험 off-target 없음 / 전체 {len(df_off)}개")

                                show_cols = ["Protein Name", "Target", "Identity(%)",
                                             "Mismatch", "E-value", "BitScore", "위험도", "Seed 위험"]
                                st.dataframe(
                                    df_off[show_cols].sort_values("Mismatch"),
                                    use_container_width=True
                                )

                                if not danger_df.empty:
                                    with st.expander(f"위험 off-target 상세 ({db_name})", expanded=True):
                                        st.dataframe(danger_df[show_cols], use_container_width=True)

                                sirna_summary[f"{db_name} 위험"] = danger_count
                                sirna_summary[f"{db_name} 전체"] = len(df_off)
                            else:
                                st.success("off-target 없음")
                                sirna_summary[f"{db_name} 위험"] = 0
                                sirna_summary[f"{db_name} 전체"] = 0

                        except Exception as e:
                            st.error(f"{db_name} 검색 실패: {e}")

                    st.markdown("#### 종합 판정")
                    if total_danger == 0:
                        st.success(
                            f"모든 DB에서 위험 off-target 없음 "
                            f"(전체 히트 {total_offtarget}개)"
                        )
                    elif total_danger <= 2:
                        st.warning(
                            f"주의 --- 위험 off-target {total_danger}개 "
                            f"(전체 히트 {total_offtarget}개)"
                        )
                    else:
                        st.error(
                            f"위험 --- off-target {total_danger}개 발견. "
                            f"다른 siRNA 후보 검토 권장"
                        )

                    sirna_summary["총 위험 off-target"] = total_danger
                    sirna_summary["총 히트 수"]         = total_offtarget
                    summary_all.append(sirna_summary)
                    st.markdown("---")

                if len(sirnas_to_check) > 1 and summary_all:
                    st.subheader("전체 siRNA 비교 요약")
                    df_summary = pd.DataFrame(summary_all)
                    st.dataframe(df_summary, use_container_width=True)

                    best = df_summary.sort_values(
                        ["총 위험 off-target", "효율 점수"],
                        ascending=[True, False]
                    ).iloc[0]
                    st.success(
                        f"최적 후보: {best['siRNA']} "
                        f"(효율 {best['효율 점수']}점 / "
                        f"위험 off-target {best['총 위험 off-target']}개)"
                    )
                    st.download_button(
                        "전체 요약 CSV 다운로드",
                        data=df_summary.to_csv(index=False),
                        file_name="offtarget_summary.csv",
                        mime="text/csv"
                    )
    raw_bp_inputs = st.text_input(
        "Enter Band Sizes for each Lane (e.g., 500, 1200, 800)", 
        value="500, 1200, 800",
        key="gel_input_field"
    )
    
    # --- 변수 초기화 및 입력 처리 (에러 방지 핵심 구간) ---
    target_bp_list = []
    if raw_bp_inputs:
        try:
            # 쉼표로 나누고 공백 제거 후 정수 변환
            target_bp_list = [int(x.strip()) for x in raw_bp_inputs.split(',') if x.strip()]
        except ValueError:
            st.error("⚠️ 숫자 형식이 올바르지 않습니다. 숫자와 쉼표만 사용해 주세요 (예: 500, 1000).")
            target_bp_list = []

    # 2. 결과 시각화 섹션 (target_bp_list가 존재할 때만 실행)
    if target_bp_list:
        with st.spinner("Generating Detailed Gel..."):
            # 레인 개수 설정 (Ladder 포함 최소 5개 레인 생성)
            num_lanes = max(5, len(target_bp_list) + 1) 
            fig, ax = plt.subplots(figsize=(num_lanes * 0.5, 3)) 
            ax.set_facecolor('#361F00') # 사진과 같은 짙은 브라운 배경
            
            # Y축 로그 스케일 및 범위 설정
            ax.set_yscale('log')
            ax.set_ylim(40, 6000)
            ax.set_xlim(0.5, num_lanes + 0.5)

            # --- 눈금자 설정 (모든 밴드 위치에 숫자 표시) ---
            ladder_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000]
            
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
            
            ax.set_yticks(ladder_sizes) 
            ax.set_yticklabels([f"{x}" for x in ladder_sizes], color='#FFA54F', fontsize=4)
            
            # X축 라벨 설정
            lane_labels = ["L"] + [f"Lane {i+1}" for i in range(num_lanes-1)]
            ax.set_xticks(range(1, num_lanes + 1))
            ax.set_xticklabels(lane_labels, color='#FFA54F', fontweight='bold', fontsize=4)
            
            # 눈금 및 테두리 스타일
            ax.tick_params(axis='y', colors='#FFA54F', length=2, labelsize=5)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.grid(False)

            # --- 메커니즘 1: 상단 Well (시료 주입구) ---
            for x in range(1, num_lanes + 1):
                well = Rectangle((x - 0.35, 4800), 0.7, 800, color='#D2691E', alpha=0.7)
                ax.add_patch(well)

            # --- 메커니즘 2: Lane L (Ladder) 밴드 그리기 ---
            for l_size in ladder_sizes:
                # 주요 지점(500, 1000, 3000)은 더 두껍게
                lw = 2.5 if l_size in [500, 1000, 3000] else 1.2
                ax.hlines(y=l_size, xmin=0.6, xmax=1.4, colors='#FF7F24', linewidth=lw, alpha=0.7)

            # --- 메커니즘 3: 각 Lane에 입력된 밴드 배정 ---
            for i, bp in enumerate(target_bp_list):
                lane_idx = i + 2 # Ladder가 1번이므로 2번 레인부터 시작
                if lane_idx <= num_lanes:
                    # 메인 밴드 (강한 주황 형광)
                    ax.hlines(y=bp, xmin=lane_idx - 0.35, xmax=lane_idx + 0.35, 
                              colors='#FF7F24', linewidth=2, alpha=0.95)
                    # 번짐 효과 (Glow)
                    ax.hlines(y=bp, xmin=lane_idx - 0.4, xmax=lane_idx + 0.4, 
                              colors='#FFA54F', linewidth=4, alpha=0.2)
                    # 밴드 위에 수치 텍스트 표시
                    ax.text(lane_idx, bp * 1.05, f"{bp}bp", color='#FF7F24', 
                            fontsize=4, ha='center', fontweight='bold', va='bottom')

            # --- 결과 출력 ---
            col_left, col_center, col_right = st.columns([1, 2, 1])
            with col_center:
                st.pyplot(fig, use_container_width=True)

            # 다운로드 버튼
            buf = io.BytesIO()
            fig.savefig(buf, format="png", facecolor='#361F00')
            st.download_button(
                label="📥 Download Realistic Gel Image", 
                data=buf.getvalue(), 
                file_name="realistic_gel_output.png", 
                mime="image/png"
            )
    else:
        st.warning("위 입력창에 bp 값을 입력해 주세요 (예: 500, 1200).")
with tab4:
    st.header("Concentration Calculator")
    st.subheader("cDNA Synthesis Calculation")
    st.info("RNA 2,000ng(2μg)을 기준으로 cDNA 합성에 필요한 RNA와 Water의 양을 계산합니다.")

    # 1. 입력 섹션 (들여쓰기 없이 깔끔하게 배치)
    rna_conc = st.number_input(
        "현재 측정된 RNA 농도를 입력하세요 (ng/μL)", 
        min_value=0.0, 
        value=481.095, 
        format="%.3f",
        help="NanoDrop 등으로 측정된 값을 입력하세요."
    )

    # 2. 계산 로직
    if rna_conc > 0:
        # RNA 볼륨 = 2000 / RNA 농도
        rna_volume = 2000 / rna_conc
        
        # Water 볼륨 = 11 - RNA 볼륨
        water_volume = 11 - rna_volume

        # 3. 결과 출력
        st.markdown("---")
        res_col1, res_col2 = st.columns(2)
        
        with res_col1:
            st.metric(label="Required RNA Volume", value=f"{rna_volume:.2f} μL")
            
        with res_col2:
            if water_volume >= 0:
                st.metric(label="Nuclease-Free Water", value=f"{water_volume:.2f} μL")
            else:
                st.error("RNA 농도가 너무 낮아 11μL를 초과합니다!")
                st.warning("RNA 양을 줄이거나 샘플을 농축해야 합니다.")

        # 4. 프로토콜 가이드 (도움말)
        if water_volume >= 0:
            st.success(f"**Total Mixture (11μL) = RNA {rna_volume:.2f}μL + Water {water_volume:.2f}μL**")
            with st.expander("상세 믹스 구성 보기"):
                st.write(f"""
                1. RNA Sample: {rna_volume:.2f} μL (2,000ng 기준)
                2. N.F. Water: {water_volume:.2f} μL
                3. **Total (Template): 11.00 μL**
                ---
                *이후 RT Master Mix 등을 추가하여 최종 볼륨을 맞추세요.*
                """)
    else:
        st.warning("농도를 입력해 주세요.")



with tab5:
    st.header("Gene 조정 및 서열 변환")
    st.write("전사, 번역, 역전사, 역번역 및 다양한 서열 변환을 수행합니다.")
    st.markdown("---")
 
    # ── 입력 섹션 ────────────────────────────────
    input_seq = st.text_area(
        "분석할 서열을 입력하세요",
        height=150,
        placeholder="DNA: ATGCGT... / RNA: AUGCGU... / Protein: MAST...",
        key="tab6_input"
    )
 
    # 입력 타입 선택
    col_type, col_codon = st.columns([1, 1])
    with col_type:
        seq_type = st.radio(
            "입력 서열 타입",
            ["DNA", "RNA", "Protein"],
            horizontal=True
        )
    with col_codon:
        codon_table = st.selectbox(
            "코돈 테이블",
            [
                "1 - Standard",
                "2 - Vertebrate Mitochondrial",
                "5 - Invertebrate Mitochondrial",
                "11 - Bacterial"
            ]
        )
        table_num = int(codon_table.split(" ")[0])
 
    st.markdown("---")
 
    # ── 변환 버튼 ────────────────────────────────
    st.subheader("변환 기능 선택")
 
    b1, b2, b3, b4, b5, b6 = st.columns(6)
    with b1:
        btn_transcribe   = st.button(" 전사\nDNA→RNA",         use_container_width=True)
    with b2:
        btn_translate    = st.button(" 번역\nDNA→Protein",     use_container_width=True)
    with b3:
        btn_rev_trans    = st.button(" 역전사\nRNA→DNA",       use_container_width=True)
    with b4:
        btn_back_trans   = st.button(" 역번역\nProtein→DNA",   use_container_width=True)
    with b5:
        btn_rev_comp     = st.button(" 역상보\nRev Complement", use_container_width=True)
    with b6:
        btn_all          = st.button(" 전체 변환",              use_container_width=True)
 
    st.markdown("---")
 
    # ── 변환 로직 ────────────────────────────────
    if input_seq.strip():
        clean_seq = "".join(input_seq.split()).upper()
 
        # 전사 (DNA → RNA)
        if btn_transcribe or btn_all:
            st.subheader(" 전사 결과 (DNA → RNA)")
            try:
                if seq_type != "DNA":
                    st.warning("전사는 DNA 입력이 필요합니다.")
                else:
                    dna    = Seq(clean_seq)
                    result = dna.transcribe()
                    st.success(f"길이: {len(result)} nt")
                    st.code(str(result), language="text")
                    st.download_button(
                        " RNA 서열 다운로드",
                        data=str(result),
                        file_name="transcription_rna.txt",
                        key="dl_transcribe"
                    )
            except Exception as e:
                st.error(f"전사 실패: {e}")
 
        # 번역 (DNA → Protein)
        if btn_translate or btn_all:
            st.subheader(" 번역 결과 (DNA → Protein)")
            try:
                if seq_type != "DNA":
                    st.warning("번역은 DNA 입력이 필요합니다.")
                else:
                    dna    = Seq(clean_seq)
                    result = dna.translate(table=table_num, to_stop=True)
                    st.success(f"길이: {len(result)} aa")
                    st.code(str(result), language="text")
 
                    # 정지코돈 포함 버전도 표시
                    result_full = dna.translate(table=table_num)
                    with st.expander("정지코돈(*) 포함 버전 보기"):
                        st.code(str(result_full), language="text")
 
                    st.download_button(
                        " 단백질 서열 다운로드",
                        data=str(result),
                        file_name="translation_protein.txt",
                        key="dl_translate"
                    )
            except Exception as e:
                st.error(f"번역 실패: 서열 길이를 확인하세요 (3의 배수). {e}")
 
        # 역전사 (RNA → DNA)
        if btn_rev_trans or btn_all:
            st.subheader(" 역전사 결과 (RNA → DNA)")
            try:
                if seq_type == "RNA":
                    rna    = Seq(clean_seq)
                    result = rna.back_transcribe()
                elif seq_type == "DNA":
                    # DNA를 RNA로 변환 후 역전사
                    rna    = Seq(clean_seq).transcribe()
                    result = rna.back_transcribe()
                    st.info("DNA 입력 → RNA 변환 후 역전사 수행")
                else:
                    st.warning("역전사는 RNA 또는 DNA 입력이 필요합니다.")
                    result = None
 
                if result:
                    st.success(f"길이: {len(result)} nt")
                    st.code(str(result), language="text")
                    st.download_button(
                        " cDNA 서열 다운로드",
                        data=str(result),
                        file_name="back_transcription_cdna.txt",
                        key="dl_rev_trans"
                    )
            except Exception as e:
                st.error(f"역전사 실패: {e}")
 
        # 역번역 (Protein → DNA 추정)
        if btn_back_trans or btn_all:
            st.subheader(" 역번역 결과 (Protein → 가능한 DNA 코돈)")
            try:
                if seq_type != "Protein":
                    st.warning("역번역은 Protein 입력이 필요합니다.")
                else:
                    # BioPython 역번역 (가장 빈도 높은 코돈 사용)
                    from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
 
                    # 각 아미노산을 대표 코돈으로 변환 (Standard table 기준)
                    codon_map = {
                        'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT',
                        'C': 'TGT', 'Q': 'CAA', 'E': 'GAA', 'G': 'GGT',
                        'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA',
                        'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT',
                        'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
                        '*': 'TAA'
                    }
                    dna_result = "".join(codon_map.get(aa, 'NNN') for aa in clean_seq)
                    result     = Seq(dna_result)
 
                    st.success(f"길이: {len(result)} nt ({len(clean_seq)} aa × 3)")
                    st.info(" 대표 코돈(most common codon) 기준으로 역번역됩니다. 실제 유전자 서열과 다를 수 있습니다.")
                    st.code(str(result), language="text")
                    st.download_button(
                        " 역번역 DNA 다운로드",
                        data=str(result),
                        file_name="back_translation_dna.txt",
                        key="dl_back_trans"
                    )
            except Exception as e:
                st.error(f"역번역 실패: {e}")
 
        # 역상보 (Reverse Complement)
        if btn_rev_comp or btn_all:
            st.subheader(" 역상보 서열 (Reverse Complement)")
            try:
                if seq_type == "Protein":
                    st.warning("역상보는 DNA 또는 RNA 입력이 필요합니다.")
                else:
                    seq    = Seq(clean_seq)
                    result = seq.reverse_complement()
                    st.success(f"길이: {len(result)} nt")
 
                    col_r1, col_r2 = st.columns(2)
                    with col_r1:
                        st.markdown("**원본 서열 (5'→3')**")
                        st.code(str(seq), language="text")
                    with col_r2:
                        st.markdown("**역상보 서열 (5'→3')**")
                        st.code(str(result), language="text")
 
                    st.download_button(
                        " 역상보 서열 다운로드",
                        data=str(result),
                        file_name="reverse_complement.txt",
                        key="dl_rev_comp"
                    )
            except Exception as e:
                st.error(f"역상보 실패: {e}")
 
        # ── 서열 기본 정보 (항상 표시) ───────────
        st.markdown("---")
        st.subheader(" 입력 서열 기본 정보")
        info_c1, info_c2, info_c3, info_c4 = st.columns(4)
 
        with info_c1:
            st.metric("서열 길이", f"{len(clean_seq)} {'aa' if seq_type == 'Protein' else 'nt'}")
        with info_c2:
            if seq_type != "Protein":
                gc = (clean_seq.count('G') + clean_seq.count('C')) / len(clean_seq) * 100 if clean_seq else 0
                st.metric("GC 함량", f"{gc:.1f}%")
            else:
                st.metric("아미노산 수", f"{len(clean_seq)} aa")
        with info_c3:
            if seq_type != "Protein":
                st.metric("A/T/G/C",
                    f"{clean_seq.count('A')}/{clean_seq.count('T')}/{clean_seq.count('G')}/{clean_seq.count('C')}"
                )
            else:
                stop_count = clean_seq.count('*')
                st.metric("정지코돈(*)", stop_count)
        with info_c4:
            if seq_type == "DNA":
                orf_count = clean_seq.count('ATG')
                st.metric("ATG 개수", orf_count)
            elif seq_type == "RNA":
                orf_count = clean_seq.count('AUG')
                st.metric("AUG 개수", orf_count)
            else:
                st.metric("서열 타입", seq_type)
 
    else:
        st.info("위 입력창에 서열을 입력하고 변환 버튼을 눌러주세요.")
     





