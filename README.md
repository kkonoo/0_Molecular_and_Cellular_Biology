# 0. Molecular and Cellular Biology

https://kkonoo.github.io/0_Molecular_and_Cellular_Biology/

HaLab에 들어온 CS·통계 전공 대학원생을 위한 분자세포생물학 Quarto 사이트입니다. 생물학 배경이 거의 없다고 가정하고, 시퀀싱 데이터(bulk·single-cell·spatial RNA-seq, ATAC-seq, variant)를 이해하는 데 필요한 만큼만 다룹니다. 구성과 배포 방식은 [1_Basic_Statistics](https://github.com/kkonoo/1_Basic_Statistics)와 같습니다.

## 파일 구조

```
.
├── _quarto.yml                  # 사이트 설정, 사이드바 목차
├── index.qmd                    # 홈 (강의 안내, 장 목록)
├── chapters/                    # 장별 .qmd (00–12, 부록 A1–A4)
├── assets/                      # CSS, AI 위젯 초기화
├── _extensions/ai-chat/         # AI 학습 도우미 위젯
└── .github/workflows/deploy.yml # main에 push하면 GitHub Pages로 배포
```

## 작성 원칙

- 교과서(Essential Cell Biology, OpenStax Biology 2e 등)는 목차 구조만 참고하고 문장·그림은 새로 씁니다.
- 그림은 직접 만든 SVG(인라인 SVG, mermaid)나 OpenStax Biology 2e(CC BY 4.0) 그림을 출처와 함께 씁니다.
- 외부 교재 링크는 4장의 [Human Genetics](https://chaek.org/books/human-genetics) 하나만 둡니다. 면역세포는 [HaLab Wiki](https://halab-wiki.knu.ac.kr/disease/immunology/index.html)로 연결합니다.
- 숫자와 명칭은 확인한 출처를 `[출처](url){.src}` 배지로 붙입니다.
- 장마다 "데이터에서는 이렇게 보여요" 상자(`::: {.callout-caution .data-box icon=false}`)를 장 끝(정리 앞)에 하나 이상 둡니다. 따로 절 제목을 달지 않습니다.
- important와 데이터 상자는 펼쳐 두고, warning·tip·note는 `collapse="true"`로 접어 둡니다.
- 코드는 데이터 상자에서 꼭 필요할 때만 짧게, base R로 씁니다. 그림에 상호작용이 필요하면 그 장에만 작은 JS를 넣습니다.
- 각 장 YAML의 `ai-context`에 그 장의 핵심 내용을 적습니다.

## 로컬 미리보기

```bash
quarto preview
```

## 처음 배포할 때

GitHub 저장소 → Settings → Pages → Source를 **GitHub Actions**로 바꿉니다.
